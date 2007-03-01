/*
 *  index_esa_base.h
 *  SeqAn
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_INDEX_ESA_BASE_H
#define SEQAN_HEADER_INDEX_ESA_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{
	
	// virtual suffix tree iterators
	template <typename TSpec = void>
	struct VSTree;


		// top down traversal iterators
		template <typename TSpec = void>
		struct TopDown;

				struct _Preorder;
				typedef Tag<_Preorder> Preorder;

			// allows an top-down iterator to go up
			template < typename TSpec = Preorder >
			struct ParentLinks {};

				// dfs order
				struct _Postorder;
				//struct _Preorder;

				typedef Tag<_Postorder> Postorder;
				//typedef Tag<_Preorder> Preorder;


		// bottom up traversal iterators
		template <typename TSpec = void>
		struct BottomUp;

			// bottom up repeat search iterators
			struct SuperMaxRepeats;
			struct SuperMaxRepeatsFast;
			struct MaxRepeats;
			struct MUMs;
			struct MaxRepeatOccurences;


/**
.Tag.Preorder:
..summary:Preorder depth-first-search.
..cat:Index
..signature:Preorder
..remarks:When given as a second parameter in @Function.goNext@ the Suffix Tree is traversed in a preorder fashion (visit the node before its children).
..see:Tag.Postorder
*/

/**
.Tag.Postorder:
..summary:Postorder depth-first-search.
..cat:Index
..signature:Postorder
..remarks:When given as a second parameter in @Function.goNext@ the Suffix Tree is traversed in a postorder fashion (visit the node after its children).
..see:Tag.Preorder
*/

/**
.Metafunction.DefaultDFSOrder:
..cat:Index
..summary:Default behaviour of @Function.goNext@ when no second parameter is given.
..signature:DefaultDFSOrder<TIterator>::Type
..param.TIterator:A @Spec.VSTree Iterator@.
..returns:$Tag.Postorder$ by default and $Tag.Preorder$ if $TIterator$ is $VSTree<TopDown<ParentLinks<> > >$ or $VSTree<TopDown<ParentLinks<Preorder> > >$.
*/
    template < typename TIterator >
    struct DefaultDFSOrder {
        typedef Postorder Type;
    };

    template < typename TIndex >
    struct DefaultDFSOrder< Iter< TIndex, VSTree< TopDown< ParentLinks<Preorder> > > > > {
        typedef Preorder Type;
    };


}

#endif
