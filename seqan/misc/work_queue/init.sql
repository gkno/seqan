CREATE TABLE jobs (
       id         INTEGER PRIMARY KEY AUTOINCREMENT,
       created_at DATETIME DEFAULT CURRENT_DATETIME,
       updated_at DATETIME DEFAULT CURRENT_DATETIME,
                  -- Point of creation and last change.
       host       VARCHAR(100),
                  -- The host the job has been started.
       status     VARCHAR(20),
                  -- The status of the job.  One of {"waiting", "running", "done"}.
       pid        INTEGER,
                  -- Process IDentifier of the job.
       tries      INTEGER DEFAULT 0,
                  -- Number of times the job has been tried to execute.
       args       TEXT
                  -- The argv array of the command.
);

CREATE TABLE job_executions (
       job_id      INTEGER REFERENCES jobs (id) ON DELETE CASCADE,
                   -- Primary key of the job this belongs to.
       created_at  DATETIME DEFAULT CURRENT_DATETIME,
                   -- Point of time of execution.
       host        VARCHAR(100),
                   -- Host the job was executed on.
       stdout      TEXT,
                   -- Stdout output.
       stderr      TEXT,
                   -- Stderr output.
       return_code INTEGER
                   -- Program return code.
);