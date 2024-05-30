Execute the next files in the order:
1. Execute Run50NeutralSimulations.sge
2. Execute RunNumberOfAllelesSimulations.sge
3. Execute RunLikelihoodInEveryFolder.sge

For example: qsub -t 1-50 Run50NeutralSimulations.sge
