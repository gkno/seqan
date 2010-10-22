#!/usr/bin/env python
import sys
import re
import os
import os.path

INVALID_IDENTIFIER = re.compile(r'\b_[A-Z_]\w*\b')
VALID_IDENTIFIERS = map(
        lambda rx: re.compile(rx),
        [ '___+',
          '^__$',
          '_N',
          '_L',
          '_H',
          '__u?int64(_t)?',
          '_FILE_OFFSET_BITS',
          '_POSIX_SYNCHRONIZED_IO',
          '__cplusplus',
          '__(force)?inline(__)?',
          '__alignof(__)?',
          '__attribute__',
          '__GLOBAL__',
          '_DELETIONS____',
          '_INSERTS______',
          '_REPLACEMENTS_',
          '__int128',
          '__SSE2__',
          '__m128i',
          '__VA_ARGS__',
          '__FILE__',
          '__LINE__',
          '__GET_OPT_H__',
          '_OPENMP',
          '__SINIX__',
          '__sgi',
          '__BEOS__',
          '__aix__',
          '__ICC',
          '__WATCOMC__',
          '__ADSPBLACKFIN__',
          '_BEOS',
          '__SUNPRO_CC?',
          '__tru64',
          '__FreeBSD__',
          '__ultrix',
          '__OPENBSD',
          '_MPRAS',
          '_HAIKU',
          '_SGI_COMPILER_VERSION',
          '_POSIX_C_SOURCE',
          '_XOPEN_SOURCE',
          '__OpenBSD__',
          '__AIX__',
          '__ADSP21000__',
          '__HAIKU__',
          '__riscos__',
          '__hpux',
          '__HP_aCC',
          '__riscos',
          '__hpua',
          '__GNUC__',
          '_ULTRIX',
          '_SCO_SV',
          '__DECCXX',
          '_XENIX',
          '__sgi__',
          '_WIN32',
          '__PGI',
          '__QNX__',
          '__APPLE__',
          '__AIX',
          '_SGI',
          '_AIX',
          '__XENIX__',
          '__INTEL_COMPILER',
          '__osf',
          '__linux__',
          '__sinix__',
          '__bsdos__',
          '__ADSPTS__',
          '__sun',
          '__sinix',
          '__NetBSD',
          '__FreeBSD',
          '__osf__',
          '__ultrix__',
          '__COMPILER_VER__',
          '__COMO__',
          '__linux',
          '__UNIX_SV__',
          '__HAIKU',
          '__WIN32__',
          '__NetBSD__',
          '__CYGWIN__',
          '_COMPILER_VERSION',
          '__BORLANDC__',
          '__TRU64__',
          '__MINGW32__',
          '__aix',
          '__BeOS',
          '__QNXNTO__',
          '__hpux__',
          '__IBMCPP__',
          '__IAR_SYSTEMS_ICC__',
          '__18CXX',
          '__HP_cc',
          '__SUNPRO_C',
          '__DECC',
          '__IBMC__',
          '_MSC_VER' ])

def valid(id):
    return any(VALID_ID.match(id) for VALID_ID in VALID_IDENTIFIERS)

def find_all(file):
    f = open(file, 'r')
    result = []
    for line in f:
        matches = INVALID_IDENTIFIER.findall(line)
        invalids = [match for match in matches if not valid(match)]
        result += invalids

    return result

def test_file_type(filename):
    pos = filename.rfind(".")
    if (pos >= 0): ext = filename[pos+1:]
    else: ext = ""

    return ext in ["c", "C", "cpp", "CPP", "c++", "C++", "h", "H", "hpp", "HPP", "h++", "H++"]

def main():
    results = {}
    project_path = sys.argv[1]
    for root, dirs, files in os.walk(project_path):
        if 'CVS' in dirs:
            dirs.remove('CVS')
        if '.svn' in dirs:
            dirs.remove('.svn')
        for file in files:
            if file.startswith('.'):
                continue # Skip hidden files.
            if file.find('_generated_forwards') != -1:
                continue # Skip generated forwards.
            path = os.path.join(root, file)
            if test_file_type(path):
                results[path] = set(find_all(path))

    all_ids = set()
    for ids in results.values():
        all_ids |= ids

    #for id in sorted(all_ids):
    #    print id
    for file in sorted(results.keys()):
        for id in results[file]:
            print '%s: %s' % (file, id)

    return 0


if __name__ == '__main__':
    sys.exit(main())
