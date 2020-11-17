def all_multiqc(wildcards):
    val = {
        'quast': all_quast(wildcards)['tsv'],
        'jellyfish': all_jellyfish(wildcards),
    }
    val['busco'] = unpack(all_busco_input(wildcards))
    print(val)
    return val
