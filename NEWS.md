# BTYD 2.4.3

 * Bug fix, thanks to Patrik Schilter. The `B1B2` component of the function `pnbd.pmf.General` defined in `R/pnbd.R` has changed. The fix is that the third parameter of the B2 function described by equation (24) [here](https://www.brucehardie.com/notes/013/Pareto_NBD_interval_pmf_rev.pdf) is now correctly set to `c <- r + s + x + 1`. The previous, incorrect version carried over the `i` component from the first term by setting `c <- w + x + 1`.

# BTYD 2.4.2

 * Now `hardie` defaults to TRUE everywhere but in `bgndb.generalParams`, where it defaults to NULL