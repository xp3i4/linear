#ifndef MAPPARM_H
#define MAPPARM_H

//efficient 
MapParm parm1 ( 
        base_block_size_,     //blockSize,
        //Const_::_DELTA,          //delta(Const_::_DELTA),
        64,                          //delta
        base_threshold_,     //threshold(Const_::_THRESHOLD),
        base_kmer_step_,       //kmerStep(Const_::_KMERSTEP),
        base_shape_len_,      //shapeLen(Const_::_SHAPELEN),
        1,                      //senstivity(0),
        0,                      //anchorDeltaThr(),
        1000,                   //minReadLen(1000),
        10,                      //listN
        20,                      //listN2
        15,                     //alpha(Const_::_ALPHA),
        5,                      //alpha2 for complex mapping 
        0.02,                    //anchorLenThr(0.02),    anchors with lenghth > this parameter is pushed into the queue
        0.5,                     //rcThr(0.75)
        0.7,                     //cordThr length of cord < cordThr are abandone
        0.7,                     //senthr: perfrom next filter on cords of length < senthr 
        0.1                      //clsthr: thread of cluster
); 

//normal
MapParm parm0 ( 
        base_block_size_,     //blockSize,
        base_delta_,          //delta(Const_::_DELTA),
        base_threshold_,     //threshold(Const_::_THRESHOLD),
        base_kmer_step_,       //kmerStep(Const_::_KMERSTEP),
        base_shape_len_,      //shapeLen(Const_::_SHAPELEN),
        0,                      //senstivity(0),
        0,                      //anchorDeltaThr(),
        1000,                   //minReadLen(1000),
        10,                      //listN
        20,                      //listN2
        15,                     //alpha(Const_::_ALPHA),
        5,                      //alpha2 for complex mapping
        0.02,                    //anchorLenThr(0.02),    anchors with lenghth > this parameter is pushed into the queue
        0.5,                     //rcThr(0.8)
        0.2,                     //cordThr length of cord < cordThr are abandoned
        0.2,                     //senthr: length of cord < senthr are erased duing path
        0.1                      //clsthr: thread of cluster

); 

//sensitive
MapParm parm2 ( 
        base_block_size_,     //blockSize,
        //Const_::_DELTA,          //delta(Const_::_DELTA),
        64,                      //delta
        base_threshold_,     //threshold(Const_::_THRESHOLD),
        base_kmer_step_,       //kmerStep(Const_::_KMERSTEP),
        base_shape_len_,      //shapeLen(Const_::_SHAPELEN),
        2,                      //senstivity(0),
        0,                      //anchorDeltaThr(),
        1000,                   //minReadLen(1000),
        4,                      //listN
        50,                      //listN2
        0.65,                     //alpha(Const_::_ALPHA),
        0.5,                      //alpha2 for complex mapping
        0.02,                    //anchorLenThr(0.02),    anchors with lenghth > this parameter is pushed into the queue
        0.8,                     //rcThr(0.75)
        0.8,                      //cordThr length of cord < cordThr are abandoned
        0.8,                     //senthr: length of cord < senthr are erased duing path
        0.1                      //clsthr: thread of cluster

        
); 


MapParm parmt ( 
        base_block_size_,     //blockSize,
        base_delta_,          //delta(Const_::_DELTA),
        base_threshold_,     //threshold(Const_::_THRESHOLD),
        base_kmer_step_,       //kmerStep(Const_::_KMERSTEP),
        base_shape_len_,      //shapeLen(Const_::_SHAPELEN),
        0,                      //senstivity(0),
        0,                      //anchorDeltaThr(),
        1000,                   //minReadLen(1000),
        2,                         //listN
        2,                      //listN2
        0.75,                     //alpha(Const_::_ALPHA),
        0.65,                      //alpha2 for complex mapping
        0.02,                    //anchorLenThr(0.02),    anchors with lenghth > this parameter is pushed into the queue
        0.8,                     //rcThr(0.8)
        0.8,                       //cordThr length of cord < cordThr are abandoned
        0.8,                     //senthr: length of cord < senthr are erased duing path
        0.1                      //clsthr: thread of cluster
        
); 

#endif
