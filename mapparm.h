// ==========================================================================
//                          Mappeing SMRT reads
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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
// Author: cxpan <chenxu.pan@fu-berlin.de>
// ==========================================================================

#ifndef MAPPARM_H
#define MAPPARM_H

//efficient 
MapParm parm1 ( 
        Const_::_BLOCKSIZE,     //blockSize,
        Const_::_DELTA,          //delta(Const_::_DELTA),
        Const_::_THRESHOLD,     //threshold(Const_::_THRESHOLD),
        Const_::_KMERSTEP,       //kmerStep(Const_::_KMERSTEP),
        Const_::_SHAPELEN,      //shapeLen(Const_::_SHAPELEN),
        //0,                      //senstivity(0),
        0,                      //anchorDeltaThr(),
        1000,                   //minReadLen(1000),
        0.75,                     //alpha(Const_::_ALPHA),
        0.65,                      //alpha2 for complex mapping 
        0.02,                    //anchorLenThr(0.02),    anchors with lenghth > this parameter is pushed into the queue
        0.5,                     //rcThr(0.75)
        0.8                     //cordThr length of cord < cordThr are abandoned
); 

//normal
MapParm parm0 ( 
        Const_::_BLOCKSIZE,     //blockSize,
        Const_::_DELTA,          //delta(Const_::_DELTA),
        Const_::_THRESHOLD,     //threshold(Const_::_THRESHOLD),
        Const_::_KMERSTEP,       //kmerStep(Const_::_KMERSTEP),
        Const_::_SHAPELEN,      //shapeLen(Const_::_SHAPELEN),
        //0,                      //senstivity(0),
        0,                      //anchorDeltaThr(),
        1000,                   //minReadLen(1000),
        0.7,                     //alpha(Const_::_ALPHA),
        0.6,                      //alpha2 for complex mapping
        0.02,                    //anchorLenThr(0.02),    anchors with lenghth > this parameter is pushed into the queue
        0.8,                     //rcThr(0.8)
        0.8                     //cordThr length of cord < cordThr are abandoned
); 

//sensitive
MapParm parm2 ( 
        Const_::_BLOCKSIZE,     //blockSize,
        Const_::_DELTA,          //delta(Const_::_DELTA),
        Const_::_THRESHOLD,     //threshold(Const_::_THRESHOLD),
        Const_::_KMERSTEP,       //kmerStep(Const_::_KMERSTEP),
        Const_::_SHAPELEN,      //shapeLen(Const_::_SHAPELEN),
        //0,                      //senstivity(0),
        0,                      //anchorDeltaThr(),
        1000,                   //minReadLen(1000),
        0.6,                     //alpha(Const_::_ALPHA),
        0.5,                      //alpha2 for complex mapping
        0.02,                    //anchorLenThr(0.02),    anchors with lenghth > this parameter is pushed into the queue
        0.8,                     //rcThr(0.75)
        0.8                      //cordThr length of cord < cordThr are abandoned
); 


MapParm parmt ( 
        Const_::_BLOCKSIZE,     //blockSize,
        Const_::_DELTA,          //delta(Const_::_DELTA),
        Const_::_THRESHOLD,     //threshold(Const_::_THRESHOLD),
        Const_::_KMERSTEP,       //kmerStep(Const_::_KMERSTEP),
        Const_::_SHAPELEN,      //shapeLen(Const_::_SHAPELEN),
        //0,                      //senstivity(0),
        0,                      //anchorDeltaThr(),
        1000,                   //minReadLen(1000),
        0.6,                     //alpha(Const_::_ALPHA),
        0.5,                      //alpha2 for complex mapping
        0.02,                    //anchorLenThr(0.02),    anchors with lenghth > this parameter is pushed into the queue
        0.8,                     //rcThr(0.8)
        0.8                       //cordThr length of cord < cordThr are abandoned
); 

#endif
