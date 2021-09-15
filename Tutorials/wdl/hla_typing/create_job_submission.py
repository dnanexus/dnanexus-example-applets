import sys
import math

input_file=sys.argv[1]
batch_size=int(sys.argv[2])

def _parse_dx_delim(delim_line):
    '''parse each list of delim output from dx find into NAME, ID, SIZE, and FOLDER'''
    id=delim_line[-1]
    split_path=delim_line[3].split('/')
    folder='/'+'/'.join(split_path[:-1])
    name=split_path[-1]
    #folder and name is not used in this example, but they can be useful for some scenerio

    return name,id,folder



fd=open(input_file)
lines=fd.readlines()
sample_number=len(lines)
batch_mapped_files=''
input_number=0
number_of_batch = int(math.ceil(sample_number*1.0/batch_size))
for batch_number in range(number_of_batch):
    batch_mapped_files=''
    for member in range(batch_size):
        delim_line = lines[input_number].strip().split('\t')
        name, id, folder = _parse_dx_delim(delim_line)
        batch_mapped_files += '-imapped_read={} '.format(id)
        final_folder='/HLA_process/' + str(batch_number)
        input_number+=1
        if input_number == sample_number:
            break

    print('dx run /arcas_hla_cram_instance_bundle \
          -ireference=genome_reference/GRCh38_full_analysis_set_plus_decoy_hla.fa {batch_mapped_files} \
          --tag 200K_exome_HLA_analysis --tag original --tag batch_n_{batch_number} \
          --folder="{final_folder}" --priority normal\
          -y --brief \
          '.format(batch_mapped_files=batch_mapped_files,batch_number=batch_number,final_folder=final_folder))
