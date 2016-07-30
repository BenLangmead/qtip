Command to make the patch:

``` cd ~/git ; make -C bwa clean ; make -C bwa-mod clean ; diff -rupN -x "Xcode" -x ".git" bwa bwa-mod > ~/git/mapq/software/bwa/bwa_conc_flags_0.7.12.patch ```