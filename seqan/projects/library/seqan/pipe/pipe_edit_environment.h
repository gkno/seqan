/*
 *  pipe_edit_environment.h
 *  SeqAn
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_PIPE_EDIT_ENVIRONMENT_H
#define SEQAN_HEADER_PIPE_EDIT_ENVIRONMENT_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

	struct _HammingDistance;
	struct _LevenshteinDistance;

	typedef Tag<_HammingDistance>		HammingDistance;
	typedef Tag<_LevenshteinDistance>	LevenshteinDistance;
	typedef Tag<_LevenshteinDistance>	EditDistance;

    template < typename TDistanceSpec >
    struct EditEnvironment;

/**
.Spec.EditEnvironment:
..cat:Pipelining
..general:Class.Pipe
..summary:Outputs tuples of the $tupleLen$ consecutive elements of the input stream.
..signature:Pipe<TInput, Tupler<tupleLen, omitLast> >
..param.TInput:The type of the pipeline module this module reads from.
..param.tupleLen:The tuple length.
...remarks:The tuples contain elements $in[i]in[i+1]...in[i+(tupleLen-1)]$.
..param.omitLast:Omit half filled tuples.
..param.omitLast:If $true$, the output stream is $tupleLen-1$ elements shorter than the input stream.
..param.omitLast:If $false$, the lengths are identical and the last tuples are filled with blanks (default constructed elements) for undefined entries.
..remarks:The output type is a @Class.Tuple@ of input elements and length $tupleLen$ (i.e. $Tuple<Value<TInput>::Type, tupleLen>$).
..remarks:The tuples are sequences of the form $in[i]in[i-1]in[i-2]..in[i-tupleLen+1]$. For $omitLast=false$ $i$ begins with 0 and for $omitLast=true$ $i$ begins with $tupleLen-1$.
*/

    //////////////////////////////////////////////////////////////////////////////
    // pipe to enumerate the hamming 1-environment
    template < typename TInput >
    struct Pipe< TInput, EditEnvironment< Tag<_HammingDistance> > >
    {
		typedef typename Value< typename Value<TInput>::Type, 2 >::Type	TTuple;
		typedef typename Value<TTuple>::Type							TValue;

        TInput                      &in;
        typename Value<Pipe>::Type	tmp, orig;
		unsigned					errorPos;		// position of substitution
		unsigned					character;		// replacement character 
		unsigned					skipChar;		// skip the original character 
        
        Pipe(TInput& _in):
            in(_in) {}

        inline typename Value<Pipe>::Type const & operator*() const {
            return tmp;
        }

        inline Pipe& operator++() {
			do {
				if (++character < ValueSize<TValue>::VALUE)
					// next replacement value
					assignValueAt(tmp.i2, errorPos, (TValue) character);
				else {
					// next error position					
					assignValueAt(tmp.i2, errorPos, orig.i2[errorPos]);
					character = 0;
					if (++errorPos < length(tmp.i2)) {
						skipChar = (unsigned) orig.i2[errorPos];
						assignValueAt(tmp.i2, errorPos, (TValue) 0);
					} else {
						// next tuple
						errorPos = 0;
						++in;
						if (!eof(in)) {
							tmp = orig = *in;
							assignValueAt(tmp.i2, 0, (TValue) 0);
						}
					}
				}
			// output the original tuple only once
			} while ((errorPos > 0) && (character == skipChar)); 

            return *this;
        }
	};


    //////////////////////////////////////////////////////////////////////////////
    // pipe to enumerate the levenshtein 1-environment
    template < typename TInput >
    struct Pipe< TInput, EditEnvironment< Tag<_LevenshteinDistance> > >
    {
		typedef typename Value< typename Value<TInput>::Type, 2 >::Type	TTuple;
		typedef typename Value<TTuple>::Type							TValue;

		enum TState { _SUBST, _DELETE, _INSERT, _INSERT_LAST, _EOF };

        TInput                      &in;
        typename Value<Pipe>::Type	tmp, orig, prev;
		unsigned					errorPos;		// position of substitution
		unsigned					character;		// replacement character 
		unsigned					skipChar;		// skip the original character 
		TState						state;

        Pipe(TInput& _in):
            in(_in),
			state(_EOF) {}

        inline typename Value<Pipe>::Type const & operator*() const {
            return tmp;
        }

        inline Pipe& operator++() {
			switch (state) {
			case _SUBST:
				// before _SUBST (tmp[1..] == orig[1..] and tmp[0] == 0) holds
				do {
					if (++character < ValueSize<TValue>::VALUE) {
						// next replacement value
						assignValueAt(tmp.i2, errorPos, (TValue) character);
					} else {
						// next substitution position
						assignValueAt(tmp.i2, errorPos, orig.i2[errorPos]);
						character = 0;
						if (++errorPos < length(tmp.i2)) {
							skipChar = (unsigned) orig.i2[errorPos];
							assignValueAt(tmp.i2, errorPos, (TValue) 0);
						} else {
							// NEXT TUPLE
							// now (tmp == orig) holds
							++in;
							if (!eos(in)) {
								prev = orig;
								orig = *in;
								tmp.i2 = orig.i2;
								assignValueAt(tmp.i2, 0, prev.i2[0]);
								assignValueAt(tmp.i2, 1, prev.i2[1]);
								if (length(tmp.i2) >= 4) {
									errorPos = 2;
									state = _DELETE;
									//::std::cout << ::std::endl << "_DELETIONS____" << ::std::endl;
									return *this;
								}
							} else {
								// LAST TUPLE
								shiftLeft(orig.i2);
								assignValueAt(tmp.i2, 0, orig.i2[0]);
								assignValueAt(tmp.i2, 1, (TValue) 0);
								character = 0;
								errorPos = 1;
								state = _INSERT_LAST;
								//::std::cout << ::std::endl << "_INSERTS______" << ::std::endl;
								return *this;
							}
						}
					}
				// output the original tuple only once
				} while ((errorPos > 0) && (character == skipChar));
				break;
			case _DELETE:
				// before _DELETE (prev=orig, ++in; tmp=orig=*in) holds
				assignValueAt(tmp.i2, errorPos, prev.i2[errorPos]);
				if (++errorPos >= length(tmp.i2) - 1) {
					assignValueAt(tmp.i2, length(tmp.i2)-1, prev.i2[length(tmp.i2)-1]);
					assignValueAt(tmp.i2, 0, orig.i2[0]);
					assignValueAt(tmp.i2, 1, (TValue) 0);
					character = 0;
					errorPos = 1;
					state = _INSERT;
					//::std::cout << ::std::endl << "_INSERTS______" << ::std::endl;
				}
				break;
			case _INSERT:
			case _INSERT_LAST:
				// before _INSERT (prev=orig, ++in; tmp=prev) holds
				if (++character < ValueSize<TValue>::VALUE)
					// next replacement value
					assignValueAt(tmp.i2, errorPos, (TValue) character);
				else {
					// next insert position					
					assignValueAt(tmp.i2, errorPos, orig.i2[errorPos]);
					character = 0;
					if (++errorPos >= length(tmp.i2) - 1 && state == _INSERT) {
						tmp = orig;
						state = _SUBST;
						//::std::cout << ::std::endl << "_REPLACEMENTS_" << ::std::endl;
						errorPos = 0;
						assignValueAt(tmp.i2, 0, (TValue) 0);
						break;
					}
					if (errorPos >= length(tmp.i2)) {
						if (eof(in))
							state = _EOF;
						else {
							tmp = orig = *in;

							// begin to insert the first char at position 0
							shiftRight(tmp.i2);
							assignValueAt(tmp.i2, 1, (TValue) 0);
							errorPos = 0;
							state = _INSERT;
							//::std::cout << ::std::endl << "_INSERTS______" << ::std::endl;
						}
						break;
					}
					assignValueAt(tmp.i2, errorPos, (TValue) 0);
				}
			default:;
			}			
            return *this;
		}
	};


    //////////////////////////////////////////////////////////////////////////////
    // global pipe functions
    template < typename TInput >
	inline bool 
	control(
		Pipe< TInput, EditEnvironment< Tag<_HammingDistance> > > &me, 
		ControlBeginRead const &command) 
	{
        if (!control(me.in, command)) return false;

		me.tmp = me.orig = *me.in;
		me.errorPos = 0;
		me.character = 0;

		return true;
	}
    
    template < typename TInput >
	inline bool 
	control(
		Pipe< TInput, EditEnvironment< Tag<_LevenshteinDistance> > > &me, 
		ControlBeginRead const &command) 
	{
        if (!control(me.in, command)) return false;

		if (eof(me.in)) {
			me.state = me._EOF;
			return true;
		}

		typedef typename Value< typename Value<TInput>::Type, 2 >::Type TTuple;
		typedef typename Value<TTuple>::Type							TValue;

		me.tmp = me.orig = *me.in;

		// begin to insert the first char at position 0
		shiftRight(me.tmp.i2);
		assignValueAt(me.tmp.i2, 0, (TValue) 0);
		me.character = 0;
		me.errorPos = 0;
		me.state = me._INSERT;
		//::std::cout << ::std::endl << "_INSERTS______" << ::std::endl;

		return true;
	}
    
    template < typename TInput >
	inline bool 
	control(
		Pipe< TInput, EditEnvironment< Tag<_LevenshteinDistance> > > &me, 
		ControlEof const &command) 
	{
		return me.state == me._EOF;
    }

    template < typename TInput >
    inline typename Size< Pipe< TInput, Pipe< TInput, EditEnvironment< Tag<_HammingDistance> > > > >::Type
    length(Pipe< TInput, EditEnvironment< Tag<_HammingDistance> > > const &me) {
		typedef typename Value< typename Value<TInput>::Type, 2 >::Type TTuple;
		unsigned alphabetSize = ValueSize< typename Value<TTuple>::Type >::VALUE;
		return length(me.in) * (1 + length(me.tmp.i2) * (alphabetSize - 1));
    }

    template < typename TInput >
    inline typename Size< Pipe< TInput, Pipe< TInput, EditEnvironment< Tag<_LevenshteinDistance> > > > >::Type
    length(Pipe< TInput, EditEnvironment< Tag<_LevenshteinDistance> > > const &me) {
		typedef typename Value< typename Value<TInput>::Type, 2 >::Type TTuple;
		unsigned alphabetSize = ValueSize< typename Value<TTuple>::Type >::VALUE;
		unsigned seqs = countSequences(me.in);
		// TODO: We run into problems when one sequence contains 1 or less tuples
		// length should be ommitted in future, but Pools or the skew algorithm needs to know the stream length
		if (length(me.in) > seqs)
			return 
				  length(me.in)      * (1 + length(me.tmp.i2) * (alphabetSize - 1)) +	// substitutions and original
				 (length(me.in) - seqs) * (length(me.tmp.i2) - 3) +						// deletions
				((length(me.in) + seqs) * (length(me.tmp.i2) - 2) + 2 * seqs) * alphabetSize;	// insertions
		else
			return 0;
    }
//}

}

#endif
