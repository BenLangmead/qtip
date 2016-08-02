Command to make the patch:

``` VER=0.7.15; cd ~/git ; make -C bwa clean ; make -C bwa-mod clean ; diff -rupN -x "Xcode" -x ".git" bwa bwa-mod > ~/git/qsim/software/bwa/bwa_conc_flags_${VER}.patch ```
