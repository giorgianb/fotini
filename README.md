# fotini

Fotini is a general-purpose, photon simulation package. It is being actively developed
for the [Figueroa Research Group](http://qit.physics.sunysb.edu/wordpress/). The purpose
of fotini is to enable the verification and debugging of optical circuits, enabling the
following workflow:
1. Design an optical circuits with the desired properties.
2. Code the optical circuit in fotini.
3. Use fotini to generate the expected output statistics.
4. Measure experimental output statistics, and see whether they are consistent with those computed by fotini.
5. Debug, using fotini as a debugging aide to identify problems in the physical setup.


Fotini works through the efficient and correct handling of [Fock States](https://en.wikipedia.org/wiki/Fock_state). 
It computes output statistics using continuous-mode quantum mechanics, and thus
for the supported optical circuit elements, reproduces all important quantum
mechanical effects. While it supports *general* spectral amplitudes for a
particular number state, it works particularly efficiently with number states
expressible as `C_1*exp(Σa_it_i^2 + Σb_i*t_i + c)`. These covers a important
variety of real-world usage.


# Example Usage
## Reproduction of Hong-Ou Mandel Effect
The Hong-Ou Mandel effect has the set-up seen in figure (1).
<img src="figs/hom-setup.png">

<figure>
  <img src="figs/hom-setup.png" alt="HOM Setup">
  <figcaption>Fig. 1 Hong-Ou Mandel Effect Setup</figcaption>
</figure>
