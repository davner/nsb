#Just a small shell wrapper to iterate over the fits files.
#Ex. nsb_process.sh kp*.fits

for file in "$@"

	do
	    ~/Downloads/code_mm/pp_prepare.py "$file"
	    ~/Downloads/code_mm/pp_register.py "$file"
	    ~/Downloads/code_mm/pp_photometry.py "$file" -aprad 6.0
	    ~/Downloads/code_mm/pp_calibrate.py "$file"
	done
