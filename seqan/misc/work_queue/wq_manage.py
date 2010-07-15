#!/usr/bin/env python
"""Work queue management.

    wq_manage init path.db
    wq_manage purge path.db
"""
from __future__ import with_statement

import optparse
import os
import os.path
import signal
import socket
import subprocess
import sys
import time

import sqlite3

USAGE_MESSAGE = """
Work Queue Management

Usage: wq_manage.py init path.db
           Initialize a new database at path.db

       wq_manage.py list path.db
           List all jobs in the database.

       wq_manage.py purge path.db
           Purge all all tasks marked as running on this machine but where no
           process with this PID is still running.

       wq_manage.py insert path.db COMMAND
           Insert the COMMAND into the database for later execution.

       wq_manage.py run path.db
           Run the jobs in the given queue.
""".strip()


def connectToDatabase(path):
    c = sqlite3.connect(path) # connect to db
    c.row_factory = dictFactory
    return c


def dictFactory(cursor, row):
    """Converts database rows into dicts, taken from Python documentation."""
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d


def initDatabase(database_path):
    if os.path.exists(database_path):
        print 'Cannot initialize database at %s.' % database_path
        print 'The file already exists!'
        return 1
    # Load init SQL.
    sql_path = os.path.join(os.path.dirname(__file__), 'init.sql')
    with open(sql_path, 'r') as f:
        sql = f.read()
    # Execute SQL.
    c = sqlite3.connect(database_path)
    cur = c.cursor()
    cur.executescript(sql)
    c.commit()
    print 'Database %s successfully initialized.' % database_path
    return 0


def pidRunning(pid):
    """Return True iff there is a process with the given pid."""
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    else:
        return True


def purgeDatabase(database_path):
    c = connectToDatabase(database_path)    # connect to db
    cur = c.cursor()
    this_hostname = socket.gethostname()  # get hostname
    # Query for the jobs running on this machine.
    sql = 'SELECT COUNT(*) FROM jobs WHERE host = ? AND status = ?'
    cur.execute(sql, (this_hostname, 'running'))
    num = cur.fetchone()['COUNT(*)']
    print 'There are %d (inexact) jobs running on this host.' % num
    sql = ('SELECT id, host, pid, status, args FROM jobs '
           'WHERE host = ? AND status = ?')
    res = cur.execute(sql, (this_hostname, 'running'))
    jobs = [job for job in res]
    # Filter to jobs not running at the moment.
    stale_jobs = [job for job in jobs if not pidRunning(job['pid'])]
    for job in stale_jobs:
        print 'Purging stale job with PID %d ("%s")...' % (job['pid'], job['args'])
        sql = 'DELETE FROM jobs WHERE id = %d'
        cur.execute(sql, (job['id'],))
    return 0


def listJobs(database_path):
    c = connectToDatabase(database_path)    # connect to db
    cur = c.cursor()
    # Query for the jobs on this machine.
    sql = 'SELECT COUNT(*) FROM jobs'
    cur.execute(sql)
    num = cur.fetchone()['COUNT(*)']
    print 'There are %d (inexact) jobs running.' % num
    sql = 'SELECT id, host, pid, status, args FROM jobs ORDER BY id'
    res = cur.execute(sql)
    print 'The following jobs are running on this machine.'
    print '%6s\t%-20s\t%-6s\t%-10s\t%s' % ('id', 'host', 'pid', 'status', 'args')
    print '-' * 70
    for job in res:
        print '%(id)6s\t%(host)-20s\t%(pid)-6s\t%(status)-10s\t%(args)s' % job
    return 0


def insertJob(database_path, command):
    c = connectToDatabase(database_path)  # connect to db
    cur = c.cursor()
    this_hostname = socket.gethostname()  # get hostname
    # Insert job into database.
    sql = ('INSERT INTO jobs (host, status, args) '
           'VALUES (?, ?, ?)')
    res = cur.execute(sql, (this_hostname, 'waiting', command))
    c.commit()
    c.close()
    return 0


def fetchJob(con):
    sql = 'SELECT COUNT(*) FROM jobs WHERE status = ?'
    res = con.execute(sql, ('waiting',))
    if int(res.fetchone()['COUNT(*)']) == 0:
        return None
    # Fetch one job.
    sql = 'SELECT id, args FROM jobs WHERE status = ? LIMIT 1'
    res = con.execute(sql, ('waiting',))
    data = res.fetchone()
    job_id, job_args = data['id'], data['args']
    # Update job's state.
    sql = 'UPDATE jobs SET status = ? WHERE id = ?'
    con.execute(sql, ('running', job_id))
    con.commit()
    return (job_id, job_args)


def executeJob(con, job_id, job_args):
    # Execute subprocess.
    print 'Running %s' % job_args
    p = subprocess.Popen(job_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # Update job entry.
    this_hostname = socket.gethostname()
    sql = 'SELECT tries FROM jobs WHERE id = ?'
    tries = con.execute(sql, (job_id,)).fetchone()['tries']
    sql = 'UPDATE jobs SET tries = ?, pid = ?, host = ? WHERE id = ?'
    con.execute(sql, (int(tries) + 1, p.pid, this_hostname, job_id))
    con.commit()
    # Wait for process to finish.
    ret = p.wait()
    # Create job_run entry.
    sql = ('INSERT INTO job_executions (job_id, host, stdout, stderr, return_code) '
           'VALUES (?, ?, ?, ?, ?)')
    con.execute(sql, (job_id, this_hostname, p.stdout.read(), p.stderr.read(), ret))
    # Update job entry.
    sql = 'UPDATE jobs SET status = ? WHERE id = ?'
    con.execute(sql, ('done', job_id))


def runJobs(database_path):
    """Fetch job after job from the database."""
    c = connectToDatabase(database_path)    # connect to db
    # Fetch all jobs.
    is_sleeping = False
    try:
        while True:
            job_id_args = fetchJob(c)
            if job_id_args:
                executeJob(c, *job_id_args)
            else:
                if not is_sleeping:
                    print 'No job. Checking in 5 s intervals.'
                    print 'Press Ctrl+C to terminate.'
                is_sleeping = True
                time.sleep(5)
            pass
    except KeyboardInterrupt:
        print 'Caught Ctrl+C. Stopping.'
        return 0


def mainInit():
    parser = optparse.OptionParser()
    options, args = parser.parse_args(sys.argv[2:])
    if len(args) != 1:
        print 'Invalid argument count!'
        print USAGE_MESSAGE
        return 1
    database_path = args[0]
    return initDatabase(database_path)


def mainPurge():
    parser = optparse.OptionParser()
    options, args = parser.parse_args(sys.argv[2:])
    if len(args) != 1:
        print 'Invalid argument count!'
        print USAGE_MESSAGE
        return 1
    database_path = args[0]
    return purgeDatabase(database_path)


def mainList():
    parser = optparse.OptionParser()
    options, args = parser.parse_args(sys.argv[2:])
    if len(args) != 1:
        print 'Invalid argument count!'
        print USAGE_MESSAGE
        return 1
    database_path = args[0]
    return listJobs(database_path)


def mainInsert():
    if len(sys.argv) < 3:
        print 'Invalid argument count!'
        print USAGE_MESSAGE
        return 1
    return insertJob(sys.argv[2], sys.argv[3:][0])


def mainRun():
    if len(sys.argv) != 3:
        print 'Invalid argument count!'
        print USAGE_MESSAGE
        return 1
    return runJobs(sys.argv[2])

def main():
    if (len(sys.argv) < 3):
        print USAGE_MESSAGE
        return 1
    command = sys.argv[1]
    if command not in ['init', 'purge', 'list', 'insert', 'run']:
        print USAGE_MESSAGE
        return 1
    FUNCS = {'init': mainInit, 'purge': mainPurge, 'list': mainList,
             'insert': mainInsert, 'run': mainRun}
    return FUNCS[command]()


if __name__ == '__main__':
    sys.exit(main())
