## gsl-splines

Wrapper functions around [`rust-GSL`]'s wrapper functions

This crates provides a cleaner interface and control over 
[`GSL`]'s splines.

Most importantly, it allows for multiple different splines to use the same Accelerator object. This provides
a significant performance boost when evaluating many splines with the same data points, at the same point, a
case which comes up a lot in many physics calculations, such as ODE problems.

[`rust-GSL`]: https://github.com/GuillaumeGomez/rust-GSL
[`GSL`]: https://www.gnu.org/software/gsl/doc/html/interp.html