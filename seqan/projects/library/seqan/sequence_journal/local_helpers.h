template< typename TSA, typename TInv >
void silly_inverse( TSA const & suffix_array, TInv & inverse_array ){
    resize( inverse_array, seqan::length(suffix_array) ); //inverse and suffix array are of same length
    for( size_t i = 0; i < length(suffix_array); ++i ){
        inverse_array[suffix_array[i]] = i; //construct inverse suffix-array in a silly way
    }
    // done . . . you happy now?
}
