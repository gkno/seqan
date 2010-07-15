CREATE TABLE jobs (
       id     INTEGER PRIMARY KEY,
       host   VARCHAR(100),
              -- The host the job has been started.
       status VARCHAR(20),
              -- The status of the job.  One of {"waiting", "running", "done"}.
       pid    INTEGER,
              -- Process IDentifier of the job.
       args   TEXT
              -- The argv array of the command.
);