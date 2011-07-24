// FRAGMENT(includes)
#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/graph_types.h>

using namespace seqan;

// FRAGMENT(tags-structs)
struct Newick_;
typedef Tag<Newick_> Newick;

struct NewickBranchLabel
{
    bool isDistanceSet;
    double distance;
    
    NewickBranchLabel() : isDistanceSet(false), distance(0)
    {}
};

// FRAGMENT(reading)
template <typename TStream, typename TSpec>
inline int read2(String<Graph<Tree<> > > & forest,
                 String<StringSet<CharString> > & vertexLabels,
                 String<String<NewickBranchLabel> > & branchLabels,
                 RecordReader<TStream, SinglePass<TSpec> > & reader,
                 Newick const & /*tag*/)
{
    return 0;
}

// FRAGMENT(writing)
template <typename TStream, typename TTree, typename TVertexDescriptor, typename TVertexLabels,
          typename TBranchLabels>
int _writeNewickRecurse(TStream & stream, TTree & tree, TVertexDescriptor v,
                        TVertexLabels & vertexLabels, TBranchLabels & branchLabels)
{
    if (numChildren(tree, v) > 0u)
    {
        int res = streamPut(stream, '(');
        if (res != 0)
            return res;
        
        typename Iterator<TTree, OutEdgeIterator>::Type it(tree, v);
        bool first = true;
        for (; !atEnd(it); goNext(it))
        {
            if (!first)
            {
                res = streamPut(stream, ',');
                if (res != 0)
                    return res;
            }
            first = false;
            res = _writeNewickRecurse(stream, tree, targetVertex(it), vertexLabels, branchLabels);
            if (res != 0)
                return res;
        }
        
        res = streamPut(stream, ')');
        if (res != 0)
            return res;
    }
    // Write label if any, quoted if required.
    if (length(property(vertexLabels, v)) > 0u)
    {
        bool needsQuoting = false;
        CharString const & label = property(vertexLabels, v);
        typename Iterator<CharString const, Rooted>::Type it = begin(label, Rooted());
        for (; !atEnd(it); ++it)
        {
            if (isblank(*it) || *it == ',' || *it == ';' || *it == '.' ||
                *it == '\'' || *it == '[' || *it == ']' || *it == '(' ||
                *it == ')')
            {
                needsQuoting = true;
                break;
            }
        }
        if (needsQuoting)
        {
            int res = streamPut(stream, '\'');
            if (res != 0)
                return res;
            it = begin(label, Rooted());
            for (; !atEnd(it); ++it)
            {
                if (*it == '\'')
                {
                    res = streamPut(stream, "''");
                    if (res != 0)
                        return res;
                }
                else
                {
                    res = streamPut(stream, *it);
                    if (res != 0)
                        return res;
                }
            }
            res = streamPut(stream, '\'');
            if (res != 0)
                return res;
        }
        else
        {
            int res = streamPut(stream, label);
            if (res != 0)
                return res;
        }
    }
    // Write branch length if any is given.
    if (property(branchLabels, v).isDistanceSet)
    {
        int res = streamPut(stream, ':');
        if (res != 0)
            return res;
        res = streamPut(stream, property(branchLabels, v).distance);
        if (res != 0)
            return res;
    }
    return 0;
}

template <typename TStream>
inline int write2(TStream & stream,
                  Graph<Tree<> > & tree,
                  StringSet<CharString> & vertexLabels,
                  String<NewickBranchLabel> & branchLabels,
                  Newick const & /*tag*/)
{
    // Write <tree>;.
    int res = _writeNewickRecurse(stream, tree, getRoot(tree), vertexLabels, branchLabels);
    if (res != 0)
        return res;
    return streamPut(stream, ';');
}

// FRAGMENT(main)
int main(int argc, char const ** argv)
{
    // Handle arguments, open file.
    if (argc != 2)
    {
        std::cerr << "Incorrect argument count!" << std::endl;
        std::cerr << "USAGE: tutorial_parse_newick INPUT.txt" << std::endl;
        return 1;
    }
    std::fstream stream(argv[1], std::ios::binary | std::ios::in);
    if (!stream.good())
    {
        std::cerr << "Could not open file " << argv[1] << std::endl;
        return 1;
    }
    RecordReader<std::fstream, SinglePass<> > reader(stream);
    
    // Load forest.
    String<Graph<Tree<> > > forest;
    String<StringSet<CharString> > vertexLabels;
    String<String<NewickBranchLabel> > branchLabels;
    int res = read2(forest, vertexLabels, branchLabels, reader, Newick());
    if (res != 0)
    {
        std::cerr << "Could not read Newick file!" << std::endl;
        return res;
    }
    
    // Dump forests.
    for (unsigned i = 0; i < length(forest); ++i)
    {
        res = write2(std::cout, forest[i], vertexLabels[i], branchLabels[i], Newick());
        std::cout << ";\n";
        if (res != 0)
        {
            std::cerr << "Error writing to stdout!" << std::endl;
            return 1;
        }
    }
    
    return 0;
}