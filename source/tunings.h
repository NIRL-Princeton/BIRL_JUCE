/*
 * tunings.h
 *
 *  Created on: Dec 18, 2019
 *      Author: josnyder
 */

#ifndef TUNINGS_H_
#define TUNINGS_H_

#include <stdint.h>
namespace birl
{
#define NUM_TUNINGS 63
        
    extern float centsDeviation[12];
    extern const float tuningPresets[NUM_TUNINGS][12];
    extern uint32_t currentTuning;
    extern uint8_t keyCenter;
//    extern const char tuningNames[NUM_TUNINGS][13];


    typedef enum _TuningNames
    {
        EqualTempered = 0,
        JustOvertone,
        Kirnberger3,
        QCMeantone,
        JI_12,
        Agricola_P,
        Kora1,
        Kora2,
        Kora3,
        Kora4,
        Johnston,
        Mambuti,
        Marimba3,
        Mbira1,
        Mbira2,
        Mbira3,
        Mbira4,
        Mbira5,
        LMY_piano,
        Xylo2,
        Xylo4,
        LMY_guitar,
        Abell1,
        Aeolic,
        Pelog,
        Slendro,
        Angklung,
        Huzam,
        Arch_chrom,
        Archytas12,
        Arch_ptol,
        Arch_sept,
        Ariel3,
        Augtetj,
        Awraamoff,
        Bagpipe1,
        Bagpipe2,
        Bagpipe3,
        Balafon,
        Balafon2,
        Balafon3,
        Balafon4,
        Bellingwolde,
        Johnston_6,
        Keenan_t9,
        Lara,
        Ligon,
        Lydian_di2,
        Lydian_di2I,
        Lydian_enh2,
        Lydian_enhI,
        Malcolms,
        Met24_chrys,
        Metals,
        Ptolemy_imix,
        Riley_albion,
        Riley_rosary,
        Couperin,
        Dowland,
        BluesMarv,
        Carlos_Harm,
        Partch_Grm,
        Partch_Greek
        } TuningNames;
    
}

#endif /* TUNINGS_H_ */
