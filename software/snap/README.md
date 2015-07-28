* Here's the command I used to make the diffs from my $HOME/git
directory on my Mac (probably should have done this with a particular
version of SNAP)

``` diff -rupN -x "Xcode" -x ".git" tmp/snap/ snap/ > ~/git/mapq/software/snap/snap_features_patch.diff ```