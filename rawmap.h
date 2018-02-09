
/*
 * serial mapping
 *
template <typename TDna, typename TSpec>
void rawMapAllComplex2(typename PMCore<TDna, TSpec>::Index   & index,
            typename PMRecord<TDna>::RecSeqs      & reads,
            typename PMRecord<TDna>::RecSeqs & genomes,
            MapParm & mapParm,
            StringSet<String<uint64_t> > & cords)
{
    typedef typename PMRecord<TDna>::RecSeq Seq;
    typedef typename PMCore<TDna, TSpec>::Anchors Anchors;
    typename PMRes::HitString hit, crhit;
    String<uint64_t> crcord;
    String<unsigned> missId;
    double time=sysTime();
    float senThr = window_size /0.85;
    MapParm complexParm = mapParm;
    complexParm.alpha = complexParm.alpha2;
    std::cerr << "complex 2\n";
    std::cerr << "Raw mapping... \r"; 

    Anchors anchors(Const_::_LLTMax, AnchorBase::size);
    Seq comStr;
    resize(cords, length(reads));
    StringSet<String<int> > f2;
    createFeatures(genomes, f2);

    for (unsigned j = 0; j < length(reads); j++)
    {
        if (length(reads[j]) < mapParm.minReadLen) // skip reads length < 1000
            continue;
        anchors.init(1);
        clear(crhit);
        mnMapReadAll<TDna, TSpec>(index, reads[j], anchors, mapParm, crhit);
        pathAll(reads[j], begin(crhit), end(crhit), f2, cords[j], mapParm.cordThr, 0);
        anchors.init(1);
        _compltRvseStr(reads[j], comStr);
        clear(crhit);
        mnMapReadAll<TDna, TSpec>(index, comStr, anchors, mapParm, crhit);
        pathAll(comStr, begin(crhit), end(crhit), f2, cords[j], mapParm.cordThr, 1);            
        shrinkToFit(cords[j]);
        if (_DefaultCord.getMaxLen(cords[j]) < length(reads[j]) / senThr
             && _DefaultCord.getMaxLen(cords[j]) > 2)
        {
            appendValue(missId, j);
            clear(cords[j]);
        }
    }
   
    for (unsigned k = 0; k < length(missId); k++)
    {
        anchors.init(1);
        if (length(reads[missId[k]]) < complexParm.minReadLen) // skip reads length < minReadLen
            continue;
        clear(crhit);
        mnMapReadAll<TDna, TSpec>(index, reads[missId[k]], anchors, complexParm, crhit);
        pathAll(reads[missId[k]], begin(crhit), end(crhit), f2, cords[missId[k]], mapParm.cordThr, 0);
        anchors.init(1);
        _compltRvseStr(reads[missId[k]], comStr);
        clear(crhit);
        mnMapReadAll<TDna, TSpec>(index, comStr, anchors, complexParm, crhit);
        pathAll(comStr, begin(crhit), end(crhit), f2, cords[missId[k]], mapParm.cordThr, 1);            
        //shrinkToFit(cords[missId[k]]);
    }
    
    std::cerr << "Raw map reads:            " << std::endl;
    std::cerr << "    End raw mapping. Time[s]: " << sysTime() - time << std::endl;
}
*/
/*
 * parallel mapping 
 */
/*
template <typename TDna, typename TSpec>
int rawMapAllComplex2(typename PMCore<TDna, TSpec>::Index   & index,
            typename PMRecord<TDna>::RecSeqs      & reads,
            typename PMRecord<TDna>::RecSeqs & genomes,
            MapParm & mapParm,
            StringSet<String<uint64_t> > & cords,
            unsigned & threads)
{
    typedef typename PMRecord<TDna>::RecSeq Seq;
    typedef typename PMCore<TDna, TSpec>::Anchors Anchors;
    typename PMRes::HitString hit;
    String<uint64_t> crcord;
    String<unsigned> missId;
    double time=sysTime();
    float senThr = window_size / mapParm.senThr;
    MapParm complexParm = mapParm;
    complexParm.alpha = complexParm.alpha2;
    complexParm.listN = complexParm.listN2;
    std::cerr << "[all_omp] Raw Mapping \n"; 

    Anchors anchors(Const_::_LLTMax, AnchorBase::size);
    Seq comStr;
    typename PMRes::HitString crhit;

    StringSet<String<int> > f2;
    double time2 = sysTime();
    createFeatures(genomes, f2, threads);
    std::cerr << "init " << sysTime() - time2 << "\n";
    //return 0;
//#pragma omp declare reduction(_joinCords : StringSet<String<uint64_t> > : omp_out = _join(omp_out, omp_in)) 
    double time1 = sysTime();
    
#pragma omp parallel//reduction(_joinCords: cords)
{
    unsigned size2 = length(reads) / threads;
    unsigned ChunkSize = size2;
    Seq comStr;
    Anchors anchors(Const_::_LLTMax, AnchorBase::size);
    typename PMRes::HitString crhit;
    StringSet<String<uint64_t> >  cordsTmp;
    String<unsigned> missIdTmp;
    unsigned thd_id =  omp_get_thread_num();
    if (thd_id < length(reads) - size2 * threads)
    {
        ChunkSize = size2 + 1;
    }
    resize(cordsTmp, ChunkSize);
    unsigned c = 0;

    #pragma omp for
    for (unsigned j = 0; j < length(reads); j++)
    {
        if (length(reads[j]) >= mapParm.minReadLen )
        {
            anchors.init(1);
            clear(crhit);
            mnMapReadList<TDna, TSpec>(index, reads[j], anchors, mapParm, crhit);
            pathAll(reads[j], begin(crhit), end(crhit), f2, cordsTmp[c], mapParm.cordThr, 0);
            anchors.init(1);
            _compltRvseStr(reads[j], comStr);
            clear(crhit);
            mnMapReadList<TDna, TSpec>(index, comStr, anchors, mapParm, crhit);
            pathAll(comStr, begin(crhit), end(crhit), f2, cordsTmp[c], mapParm.cordThr, 1);            
           if (_DefaultCord.getMaxLen(cordsTmp[c]) < length(reads[j]) / senThr)// && 
               //_DefaultCord.getMaxLen(cordsTmp[c]) > 0)
           {
               appendValue(missIdTmp, j);
               clear(cordsTmp[c]);
           }   
        }
        c += 1;
    } 
    #pragma omp for ordered
    for (unsigned j = 0; j < threads; j++)
        #pragma omp ordered
        {
            append(cords, cordsTmp);
            append(missId, missIdTmp);
        }
}

    std::cerr << "map1 " << sysTime() - time1 << std::endl;
    time1 = sysTime();
    
    for (unsigned k = 0; k < length(missId); k++)
    {

//        std::cout << "[miis] " << missId[k] << std::endl;
        anchors.init(1);
        if (length(reads[missId[k]]) < complexParm.minReadLen) // skip reads length < minReadLen
            continue;
        clear(crhit);
        clear(cords[missId[k]]);
        mnMapReadList<TDna, TSpec>(index, reads[missId[k]], anchors, complexParm, crhit);
        pathAll(reads[missId[k]], begin(crhit), end(crhit), f2, cords[missId[k]], complexParm.cordThr, 0);

        anchors.init(1);
        _compltRvseStr(reads[missId[k]], comStr);
        clear(crhit);
        mnMapReadList<TDna, TSpec>(index, comStr, anchors, complexParm, crhit);
        pathAll(comStr, begin(crhit), end(crhit), f2, cords[missId[k]], complexParm.cordThr, 1);            
        //shrinkToFit(cords[missId[k]]);
    }
    
    std::cerr << "map2 " << sysTime() - time1 << std::endl;

    std::cerr << "Raw map reads:            " << length(missId) << std::endl;
    std::cerr << "    End raw mapping. Time[s]: " << sysTime() - time << std::endl;
    return 0;
}

/*
 * parall mapping, this version also parallized the part of the secound round
 * mapping. Hoever the senstivity is lower than the the last one. The reason is 
 * unchecked.
 *
template <typename TDna, typename TSpec>
int rawMapAllComplex2(typename PMCore<TDna, TSpec>::Index   & index,
            typename PMRecord<TDna>::RecSeqs      & reads,
            typename PMRecord<TDna>::RecSeqs & genomes,
            MapParm & mapParm,
            StringSet<String<uint64_t> > & cords,
            unsigned & threads)
{
    typedef typename PMRecord<TDna>::RecSeq Seq;
 
    double time=sysTime();
    float senThr = window_size / mapParm.senThr;
    MapParm complexParm = mapParm;
    complexParm.alpha = complexParm.alpha2;
    complexParm.listN = complexParm.listN2;
    std::cerr << "[all_omp] Raw Mapping \n"; 


    StringSet<String<int> > f2;
    double time2 = sysTime();
    createFeatures(genomes, f2, threads);
    std::cerr << "init " << sysTime() - time2 << "\n";
    //return 0;
//#pragma omp declare reduction(_joinCords : StringSet<String<uint64_t> > : omp_out = _join(omp_out, omp_in)) 
    double time1 = sysTime();
    
#pragma omp parallel//reduction(_joinCords: cords)
{
    unsigned size2 = length(reads) / threads;
    unsigned ChunkSize = size2;
    Seq comStr;
    Anchors anchors(Const_::_LLTMax, AnchorBase::size);
    typename PMRes::HitString crhit;
    StringSet<String<uint64_t> >  cordsTmp;
    unsigned thd_id =  omp_get_thread_num();
    if (thd_id < length(reads) - size2 * threads)
    {
        ChunkSize = size2 + 1;
    }
    resize(cordsTmp, ChunkSize);
    unsigned c = 0;

    #pragma omp for
    for (unsigned j = 0; j < length(reads); j++)
    {
        if (length(reads[j]) >= mapParm.minReadLen )
        {
            anchors.init(1);
            clear(crhit);
            mnMapReadList<TDna, TSpec>(index, reads[j], anchors, mapParm, crhit);
            pathAll(reads[j], begin(crhit), end(crhit), f2, cordsTmp[c], mapParm.cordThr, 0);
            
            anchors.init(1);
            _compltRvseStr(reads[j], comStr);
            clear(crhit);
            mnMapReadList<TDna, TSpec>(index, comStr, anchors, mapParm, crhit);
            pathAll(comStr, begin(crhit), end(crhit), f2, cordsTmp[c], mapParm.cordThr, 1);            
            if (_DefaultCord.getMaxLen(cordsTmp[c]) < length(reads[j]) / senThr)// && 
               //_DefaultCord.getMaxLen(cordsTmp[c]) > 0)
            {
                clear(cordsTmp[c]);
                
                anchors.init(1);
                clear(crhit);
                mnMapReadList<TDna, TSpec>(index, reads[j], anchors, complexParm, crhit);
                pathAll(reads[j], begin(crhit), end(crhit), f2, cordsTmp[c], complexParm.cordThr, 0);
               
                anchors.init(1);
                clear(crhit);
                mnMapReadList<TDna, TSpec>(index, comStr, anchors, complexParm, crhit);
                pathAll(comStr, begin(crhit), end(crhit), f2, cordsTmp[c], complexParm.cordThr, 1);
               
            }   
        }
        c += 1;
    } 
    #pragma omp for ordered
    for (unsigned j = 0; j < threads; j++)
        #pragma omp ordered
        {
            append(cords, cordsTmp);
        }
}

    std::cerr << "    End raw mapping. Time[s]: " << sysTime() - time << std::endl;
    return 0;
}
*/
/*
template <typename TDna, typename TSpec>
void rawMapFirstComplex2(typename PMCore<TDna, TSpec>::Index   & index,
            typename PMRecord<TDna>::RecSeqs      & reads,
            typename PMRecord<TDna>::RecSeqs & genomes,
            MapParm & mapParm,
            StringSet<String<uint64_t> > & cords,
            unsigned & threads)
{
    typedef typename PMRecord<TDna>::RecSeq Seq;
    typedef typename PMCore<TDna, TSpec>::Anchors Anchors;
    typename PMRes::HitString hit;
    String<uint64_t> crcord;
    String<unsigned> missId;
    double time=sysTime();
    float senThr = window_size /0.85;
    MapParm complexParm = mapParm;
    complexParm.alpha = complexParm.alpha2;
    std::cerr << "[first] Raw Mapping \n"; 

    Anchors anchors(Const_::_LLTMax, AnchorBase::size);
    Seq comStr;
    typename PMRes::HitString crhit;

    StringSet<String<int> > f2;
    createFeatures(genomes, f2);

    double time1 = sysTime();
#pragma omp parallel
{
    unsigned size2 = length(reads) / threads;
    unsigned ChunkSize = size2;
    Seq comStr;
    Anchors anchors(Const_::_LLTMax, AnchorBase::size);
    typename PMRes::HitString crhit;
    StringSet<String<uint64_t> >  cordsTmp;
    String<unsigned> missIdTmp;
    unsigned thd_id =  omp_get_thread_num();
    if (thd_id < length(reads) - size2 * threads)
    {
        ChunkSize = size2 + 1;
    }
    resize(cordsTmp, ChunkSize);
    unsigned c = 0;

    #pragma omp for
    for (unsigned j = 0; j < length(reads); j++)
    {
        if (length(reads[j]) >= mapParm.minReadLen )
        {
            anchors.init(1);
            clear(crhit);
            mnMapReadFirst<TDna, TSpec>(index, reads[j], anchors, mapParm, crhit);
            pathAll(reads[j], begin(crhit), end(crhit), f2, cordsTmp[c], mapParm.cordThr, 0);
            anchors.init(1);
            _compltRvseStr(reads[j], comStr);
            clear(crhit);
            mnMapReadFirst<TDna, TSpec>(index, comStr, anchors, mapParm, crhit);
            pathAll(comStr, begin(crhit), end(crhit), f2, cordsTmp[c], mapParm.cordThr, 1);            
           if (_DefaultCord.getMaxLen(cordsTmp[c]) < length(reads[j]) / senThr &&
               _DefaultCord.getMaxLen(cordsTmp[c]) >= 0)
           {
               printf("[debug] maxlen %d\n", _DefaultCord.getMaxLen(cordsTmp[c]));
               appendValue(missIdTmp, j);
               clear(cordsTmp[c]);
           }   
        }
        c += 1;
    } 
    #pragma omp for ordered
    for (unsigned j = 0; j < threads; j++)
        #pragma omp ordered
        {
            append(cords, cordsTmp);
            append(missId, missIdTmp);
        }
}

    std::cerr << "map1 " << sysTime() - time1 << std::endl;
    time1 = sysTime();
    
    for (unsigned k = 0; k < length(missId); k++)
    {

        anchors.init(1);
        if (length(reads[missId[k]]) < complexParm.minReadLen) // skip reads length < minReadLen
            continue;
        clear(crhit);
        clear(cords[missId[k]]);
        mnMapReadAll<TDna, TSpec>(index, reads[missId[k]], anchors, complexParm, crhit);
        pathAll(reads[missId[k]], begin(crhit), end(crhit), f2, cords[missId[k]], mapParm.cordThr, 0);

        anchors.init(1);
        _compltRvseStr(reads[missId[k]], comStr);
        clear(crhit);
        mnMapReadAll<TDna, TSpec>(index, comStr, anchors, complexParm, crhit);
        pathAll(comStr, begin(crhit), end(crhit), f2, cords[missId[k]], mapParm.cordThr, 1);            
        shrinkToFit(cords[missId[k]]);
    }
    
    
    std::cerr << "map2 " << sysTime() - time1 << std::endl;

    std::cerr << "Raw map reads:            " << length(missId) << std::endl;
    std::cerr << "    End raw mapping. Time[s]: " << sysTime() - time << std::endl;
} 
