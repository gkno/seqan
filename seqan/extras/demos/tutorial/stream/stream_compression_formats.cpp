// FRAGMENT(header)
#include <cstdio>
#include <fstream>
#if SEQAN_HAS_ZLIB
#include <zlib.h>
#endif  // #if SEQAN_HAS_ZLIB
#if SEQAN_HAS_BZIP2
#include <bzlib.h>
#endif  // #if SEQAN_HAS_BZIP2

#include <seqan/basic.h>
#include <seqan/stream.h>

using namespace seqan;

// FRAGMENT(open-gz)
int openGz(char const * filename)
{
#if SEQAN_HAS_ZLIB
    gzFile f = gzopen(filename, "rb");
    if (gzdirect(f))
    {
        std::cerr << "ERROR: GZip file has the wrong format!" << std::endl;
        gzclose(f);
        return 1;
    }
    Stream<GZFile> f2(f);
    
    // Totally inefficient char-wise writing of characters from .gz file to stderr.
    while (!streamEof(f2))
    {
        char c = '\0';
        int res = streamReadChar(c, f2);
        if (res != 0)
        {
            std::cerr << "ERROR: Reading byte from GZip file." << std::endl;
            return 1;
        }
        std::cout << c;
    }
    
    gzclose(f);
#else  // #if SEQAN_HAS_ZLIB
    std::cerr << "ZLIB not available!" << std::endl;
#endif  // #if SEQAN_HAS_ZLIB
    return 0;
}

// FRAGMENT(open-bz2)
int openBz2(char const * filename)
{
#if SEQAN_HAS_BZIP2
    FILE * f = fopen(filename, "rb");
    if (!f)
    {
        std::cerr << "ERROR: Could not open input file!" << std::endl;
        return 1;
    }
    int err = BZ_OK;
    BZFILE * f2 = BZ2_bzReadOpen(&err, f, 0, 0, NULL, 0);
    if (err != BZ_OK)
    {
        std::cerr << "ERROR: Could not open bzip2 file!" << std::endl;
        BZ2_bzReadClose(&err, f2);
        if (err != BZ_OK)
            std::cerr << "ERROR: Could not close bzip2 file!" << std::endl;
        fclose(f);
    }
    Stream<BZ2File> f3(f2);

    // Totally inefficient char-wise writing of characters from .bz2 file to stderr.
    while (!streamEof(f3))
    {
        char c = '\0';
        int res = streamReadChar(c, f3);
        if (res != 0)
        {
            std::cout << "ERROR: Reading byte from BZ2 file." << std::endl;
            return 1;
        }
        std::cerr << c;
    }
    
    BZ2_bzReadClose(&err, f2);
    if (err != BZ_OK)
    {
        std::cerr << "ERROR: Could not close bzip2 file!" << std::endl;
        fclose(f);
    }
    fclose(f);
#else  // #if SEQAN_HAS_BZIP2
    std::cerr << "BZLIB not available!" << std::endl;
#endif  // #if SEQAN_HAS_BZIP2
    return 0;
}

// FRAGMENT(main)
int main(int argc, char const ** argv)
{
    if (argc != 2)
        return 1;
    openGz(argv[1]);
    openBz2(argv[1]);
    return 0;
}
