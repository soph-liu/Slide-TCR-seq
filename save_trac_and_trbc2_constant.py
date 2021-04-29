import pysam
import pickle as pkl

directory_start_macosko = '/broad/macosko/'
def save_trac_and_trbc2_constant_umis(puck_name, species):
    abbrev_puck_name = '_'.join(puck_name.split('_')[1:])
    
    fn = '{}data/libraries/{}/{}.bam'.format(directory_start_macosko, puck_name, abbrev_puck_name)

    samfile = pysam.AlignmentFile(fn, "rb")
    
    if species == 'mouse':
        tcrb = 'Trbc2'
        tcra = 'Trac'
    if species == 'human':
        tcrb = 'TRBC2'
        tcra = 'TRAC'

    save_trac_trbc2 = {}
    save_trac_trbc2[tcra] = []
    save_trac_trbc2[tcrb] = []

    read_ct = 0
    for read in samfile.fetch():
        read_ct += 1
        if read_ct % 1000000 == 0:
            print(read_ct)
        if read.has_tag('XG'):
            if read.get_tag('XG') == tcra:
                cell_barcode = read.get_tag('XC')
                umi = read.get_tag('XM')
                save_trac_trbc2[tcra].append([cell_barcode,umi])
            if read.get_tag('XG') == tcrb:
                cell_barcode = read.get_tag('XC')
                umi = read.get_tag('XM')
                save_trac_trbc2[tcrb].append([cell_barcode,umi])
    
    with open('/broad/thechenlab/Fei/TCR/Wu_200909/trac_trbc2_constant_umis_{}.pickle'.format(puck_name), 'wb') as handle:
        pkl.dump(save_trac_trbc2, handle, protocol=pkl.HIGHEST_PROTOCOL)

    return save_trac_trbc2

save_trac_and_trbc2_constant_umis('2020-12-14_Puck_200528_22_s','human')
