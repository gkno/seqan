#ifndef BENCHMARKS_READ_MAPPERS_RETURN_CODES_H_
#define BENCHMARKS_READ_MAPPERS_RETURN_CODES_H_

// Define some return codes.
const int kRetOk = 0;       // OK, no errors.
const int kRetArgsErr = 1;  // Errors in arguments.
const int kRetIoErr = 2;    // I/O error, problem reading files.
const int kFatalErr = 3;    // Some other sort of fatal error.

#endif  // BENCHMARKS_READ_MAPPERS_RETURN_CODES_H_
