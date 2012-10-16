#!/usr/bin/env python

import subprocess
import dxpy

config_template = '''
<Files>
{Filenames}

<Options>
ThreadNumber={ThreadNumber}
PairedEnd={PairedEnd} // Possible values: True, False. Default value=False
FileFormat={FileFormat} // Possible values: FASTQ, FASTA, QSEQ. Default value=FASTQ
AutoPenalty={AutoPenalty} // Possible values: True, False. Default value=True
FixedPenalty={FixedPenalty} // Possible values: 0-100. Default value=2, when AutoPenalty=False
Gzip={Gzip} // Possible values: True, False. Default value=False
ExpressionMeasurement={ExpressionMeasurement} // possible values: None, TPM, RPKM, TPM_Transcript, RPKM_Transcript. Default value =None
SearchNovelExonJunction={SearchNovelExonJunction} // Possible values: True, False. Default value=False
GenerateSamFiles={GenerateSamFiles} // Possible values: True, False. Default value=False

<Output>
OutputName=output
OutputPath=.
'''

@dxpy.entry_point('main')
def main(**kwargs):
    filenames = []
    for fastq_file in kwargs['fastq_files']:
        name = dxpy.describe(fastq_file)['name']
        dxpy.download_dxfile(fastq_file, name)
        name = fastq_file
        filenames.append("./" + name)
    with open("secontrol_linux.ini", "w") as fh:
        fh.write(config_template.format(Filenames="\n".join(filenames),
                                        ThreadNumber=kwargs.get("ThreadNumber", 4),
                                        PairedEnd=kwargs.get("PairedEnd", False),
                                        FileFormat=kwargs.get("FileFormat", "FASTQ"),
                                        AutoPenalty=kwargs.get("AutoPenalty", True),
                                        FixedPenalty=kwargs.get("FixedPenalty", 2),
                                        Gzip=kwargs.get("Gzip", False),
                                        ExpressionMeasurement=kwargs.get("ExpressionMeasurement", None),
                                        SearchNovelExonJunction=kwargs.get("SearchNovelExonJunction", True),
                                        GenerateSamFiles=kwargs.get("GenerateSamFiles", True)))
    subprocess.call("mono ./resources/opt/OSAv1.8.1/osa.exe --alignrna ./resources/opt/OSAv1.8.1 Human.B37 RefGene secontrol_linux.ini",
                    shell=True)
