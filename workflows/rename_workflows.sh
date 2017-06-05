#!/bin/bash
find . -type f -print0 | xargs -0 sed -i 's/milkyway_msgfplus_search/milkyway_msgfplus_search/g'
find . -type f -print0 | xargs -0 sed -i 's/milkyway_crux_percolator/milkyway_crux_percolator/g'
find . -type f -print0 | xargs -0 sed -i 's/milkyway_percolator_to_fido/milkyway_percolator_to_fido/g'
find . -type f -print0 | xargs -0 sed -i 's/milkyway_crux_pin_fixer/milkyway_crux_pin_fixer/g'
find . -type f -print0 | xargs -0 sed -i 's/milkyway_msgf2pin_converter/milkyway_msgf2pin_converter/g'
find . -type f -print0 | xargs -0 sed -i 's/milkyway_spectral_counts/milkyway_spectral_counts/g'
