// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// This is the umbrella header for the graph_decomposition module.
//
// This module contains algorithms for clustering and decomposing graphs.
// Graph decomsition can be seen as clustering with more structural
// information.
// ==========================================================================

#ifndef SEQAN_SCORE_H_
#define SEQAN_SCORE_H_

// ==========================================================================
// Prerequisites
// ==========================================================================

#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
#include <seqan/misc/misc_bit_twiddling.h>

// ==========================================================================
// Generated Forwards
// ==========================================================================

#ifdef SEQAN_SWITCH_USE_FORWARDS
//#include <seqan/graph_decomposition/graph_decomposition_generated_forwards.h>
#endif

// ==========================================================================
// Module Headers
// ==========================================================================

// --------------------------------------------------------------------------
// Basic Definitions
// --------------------------------------------------------------------------

#include <seqan/graph_decomposition/graph_decomposition_base.h>

// --------------------------------------------------------------------------
// Clustering Algorithms
// --------------------------------------------------------------------------

#include <seqan/graph_decomposition/graph_cluster_digraph_tarjan.h>

// --------------------------------------------------------------------------
// Decomposition Algorithms
// --------------------------------------------------------------------------

#include <seqan/graph_decomposition/graph_decompose_graph_stiege.h>
#include <seqan/graph_decomposition/graph_decompose_digraph_condensed.h>
#include <seqan/graph_decomposition/graph_decompose_digraph_stiege.h>

#endif  // SEQAN_SCORE_H_
