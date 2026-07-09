# make dir
mkdir -p ./data
mkdir -p ./figure
mkdir -p ./bin

python gap_eval.py plot
python check_calibration.py

