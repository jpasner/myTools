mkdir output_histos_from_container
for container in $(ls /afs/cern.ch/work/j/jpasner/public/HbbISR/samples/)
do
  for sample in $(ls /afs/cern.ch/work/j/jpasner/public/HbbISR/samples/$container)
  do
    mkdir output_histos_from_container/$container
    xAH_run.py --files /afs/cern.ch/work/j/jpasner/public/HbbISR/samples/$container/$sample --submitDir output_histos_from_container/$container/$sample --config ../HbbISR/data/config_fatjet.py --treeName outTree direct
  done
done
