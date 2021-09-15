version 1.0

task arcas_hla_cram_instance_bundle{
    input {
        Array[File]+ mapped_read
        File reference
    }

    command <<<
        set -x -e -o pipefail
        mkdir output
        for input in ~{sep=" " mapped_read}; do
            file_prefix=$( basename $input ".cram")
            time samtools view -b -h ${input} -T ~{reference} > ${file_prefix}.bam
            time arcasHLA extract ${file_prefix}.bam -o output --paired -t $(nproc) -v
            time arcasHLA genotype output/${file_prefix}.extracted.1.fq.gz output/${file_prefix}.extracted.2.fq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o output -t 8 -v
            rm output/*.fq.gz output/*.alignment.p ${file_prefix}.bam
        done
    >>>
    output {
        Array[File] genotype = glob("output/*.genotype.json")
        Array[File] gene = glob("output/*.genes.json")
        Array[File] log_files = glob("output/*.log")
    }
    runtime {
        docker: "dx://hla_project:/Docker/arcas_hla_0.0.1.tar.gz"
        dx_timeout: "48H"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    parameter_meta {
    mapped_read: {
        description: "mapped short read data",
        patterns: ["*.cram"],
        stream: true
    }
    reference: {
        description: "reference genome",
        patterns: ["*.fa","*.fasta"]
    }
    }
}
