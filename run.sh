set -e
python3 -m experiments.run_asymptotics
python3 -m experiments.save_tables
python3 -m experiments.make_plots
echo "Artifacts in ./figures"
