.. sirocco documentation master file, created by
   sphinx-quickstart on Sun Jan 14 18:04:35 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

####################################################################################
SIROCCO - Simulating Ionization and Radiation in Outflows Created by Compact Objects
####################################################################################

figure:: ../images/logo.png

SIROCCO is a Monte-Carlo radiative transfer code designed to simulate the spectrum of biconical (or spherical)
winds in disk systems.  It was formerly known as Python, and originally written by
`Long and Knigge (2002) <https://ui.adsabs.harvard.edu/abs/2002ApJ...579..725L/abstract>`_ and
was intended for simulating the spectra of winds in cataclysmic variables. Since then, it has
also been used to simulate the spectra of systems ranging from young stellar objects to AGN. 
SIROCCO is named after the `Sirocco wind <https://en.wikipedia.org/wiki/Sirocco>`_, and also 
stands for Simulating Ionization and Radiation in Outflows Created by Compact Objects. 

The program is written in C and can be compiled on systems runining various flavors of linux, including macOS and the
Windows Subsystem for Linux (WSL). The code is is available on `GitHub <https://github.com/sirocco-rt/sirocco>`_. Issues
regarding the code and suggestions for improvement should be reported there.  We actively encourage other to make use of
the code for their own science.  If anyone has questions about whether the code might be useful for a project, we
encourage you to contact one of the authors of the code.

---------------------------------------
Documentation \& Publications
---------------------------------------

Various documentation exists:

* A :doc:`Quick Guide <quick>` describing how to install and run SIROCCO (in a fairly mechanistic fashion).
* More detailed documentation on this site
* A code release paper, submitted in October 2024
* Various PhD theses that describe 

For more information on how this page was generated and how to create documentation for SIROCCO,
look at the page for :doc:`documentation on the documentation <meta>`.

The following papers have used SIROCCO: 

- Parkinson, E. J., Knigge, C., Dai, L., Thomsen, L. L., Matthews, J. H., & Long, K. S. (2024). *arXiv*, arXiv:2408.16371. DOI: [10.48550/arXiv.2408.16371](https://doi.org/10.48550/arXiv.2408.16371), ADS: [Link](https://ui.adsabs.harvard.edu/abs/2024arXiv240816371P)

- Tampo, Y., Knigge, C., Long, K. S., Matthews, J. H., & Segura, N. C. (2024). *MNRAS*, 532, 1199. DOI: [10.1093/mnras/stae1557](https://doi.org/10.1093/mnras/stae1557), ADS: [Link](https://ui.adsabs.harvard.edu/abs/2024MNRAS.532.1199T)

- Higginbottom, N., Scepi, N., Knigge, C., Long, K. S., Matthews, J. H., & Sim, S. A. (2024). *MNRAS*, 527, 9236. DOI: [10.1093/mnras/stad3830](https://doi.org/10.1093/mnras/stad3830), ADS: [Link](https://ui.adsabs.harvard.edu/abs/2024MNRAS.527.9236H)

- Matthews, J. H., Strong-Wright, J., Knigge, C., Hewett, P., Temple, M. J., Long, K. S., Rankine, A. L., et al. (2023). *MNRAS*, 526, 3967. DOI: [10.1093/mnras/stad2895](https://doi.org/10.1093/mnras/stad2895), ADS: [Link](https://ui.adsabs.harvard.edu/abs/2023MNRAS.526.3967M)

- Koljonen, K. I. I., Long, K. S., Matthews, J. H., & Knigge, C. (2023). *MNRAS*, 521, 4190. DOI: [10.1093/mnras/stad809](https://doi.org/10.1093/mnras/stad809), ADS: [Link](https://ui.adsabs.harvard.edu/abs/2023MNRAS.521.4190K)

- Parkinson, E. J., Knigge, C., Matthews, J. H., Long, K. S., Higginbottom, N., & Sim, S. A. (2022). *MNRAS*, 510, 5426. DOI: [10.1093/mnras/stac027](https://doi.org/10.1093/mnras/stac027), ADS: [Link](https://ui.adsabs.harvard.edu/abs/2022MNRAS.510.5426P)

- Middleton, M. J., Higginbottom, N., Knigge, C., Khan, N., & Wiktorowicz, G. (2022). *MNRAS*, 509, 1119. DOI: [10.1093/mnras/stab2991](https://doi.org/10.1093/mnras/stab2991), ADS: [Link](https://ui.adsabs.harvard.edu/abs/2022MNRAS.509.1119M)

- Parkinson, E. J., Knigge, C., Long, K. S., Matthews, J. H., Higginbottom, N., Sim, S. A., & Hewitt, H. A. (2020). *MNRAS*, 494, 4914. DOI: [10.1093/mnras/staa1060](https://doi.org/10.1093/mnras/staa1060), ADS: [Link](https://ui.adsabs.harvard.edu/abs/2020MNRAS.494.4914P)

- Matthews, J. H., Knigge, C., Higginbottom, N., Long, K. S., Sim, S. A., Mangham, S. W., Parkinson, E. J., et al. (2020). *MNRAS*, 492, 5540. DOI: [10.1093/mnras/staa136](https://doi.org/10.1093/mnras/staa136), ADS: [Link](https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.5540M)

- Higginbottom, N., Knigge, C., Sim, S. A., Long, K. S., Matthews, J. H., & Hewitt, H. A. (2020). *MNRAS*, 492, 5271. DOI: [10.1093/mnras/staa209](https://doi.org/10.1093/mnras/staa209), ADS: [Link](https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.5271H)

- Charles, P., Matthews, J. H., Buckley, D. A. H., Gandhi, P., Kotze, E., & Paice, J. (2019). *MNRAS*, 489, L47. DOI: [10.1093/mnrasl/slz120](https://doi.org/10.1093/mnrasl/slz120), ADS: [Link](https://ui.adsabs.harvard.edu/abs/2019MNRAS.489L..47C)

- Mangham, S. W., Knigge, C., Williams, P., Horne, K., Pancoast, A., Matthews, J. H., Long, K. S., et al. (2019). *MNRAS*, 488, 2780. DOI: [10.1093/mnras/stz1713](https://doi.org/10.1093/mnras/stz1713), ADS: [Link](https://ui.adsabs.harvard.edu/abs/2019MNRAS.488.2780M)

- Higginbottom, N., Knigge, C., Long, K. S., Matthews, J. H., & Parkinson, E. J. (2019). *MNRAS*, 484, 4635. DOI: [10.1093/mnras/stz310](https://doi.org/10.1093/mnras/stz310), ADS: [Link](https://ui.adsabs.harvard.edu/abs/2019MNRAS.484.4635H)

- Milliner, K., Matthews, J. H., Long, K. S., & Hartmann, L. (2019). *MNRAS*, 483, 1663. DOI: [10.1093/mnras/sty3197](https://doi.org/10.1093/mnras/sty3197), ADS: [Link](https://ui.adsabs.harvard.edu/abs/2019MNRAS.483.1663M)

- Higginbottom, N., Knigge, C., Long, K. S., Matthews, J. H., Sim, S. A., & Hewitt, H. A. (2018). *MNRAS*, 479, 3651. DOI: [10.1093/mnras/sty1599](https://doi.org/10.1093/mnras/sty1599), ADS: [Link](https://ui.adsabs.harvard.edu/abs/2018MNRAS.479.3651H)

- Mangham, S. W., Knigge, C., Matthews, J. H., Long, K. S., Sim, S. A., & Higginbottom, N. (2017). *MNRAS*, 471, 4788. DOI: [10.1093/mnras/stx1863](https://doi.org/10.1093/mnras/stx1863), ADS: [Link](https://ui.adsabs.harvard.edu/abs/2017MNRAS.471.4788M)

- Higginbottom, N., Proga, D., Knigge, C., & Long, K. S. (2017). *ApJ*, 836, 42. DOI: [10.3847/1538-4357/836/1/42](https://doi.org/10.3847/1538-4357/836/1/42), ADS: [Link](https://ui.adsabs.harvard.edu/abs/2017ApJ...836...42H)

- Matthews, J. H. (2016). *PhDT*. DOI: [10.5281/zenodo.1256805](https://doi.org/10.5281/zenodo.1256805), ADS: [Link](https://ui.adsabs.harvard.edu/abs/2016PhDT.......348M)

- Matthews, J. H., Knigge, C., Long, K. S., Sim, S. A., Higginbottom, N., & Mangham, S


-------
Authors
-------
The authors of the SIROCCO code and their institutions are:

Knox Long
  Space Telescope Science Institute, 3700 San Martin Drive, Baltimore, MD 21218, USA
  Eureka Scientific, Inc., 2452 Delmer St., Suite 100, Oakland, CA 94602-3017, USA

Christian Knigge
  Department of Physics and Astronomy, University of Southampton, Southampton, SO17 1BJ, UK

Stuart Sim
  School of Mathematics and Physics, Queen's University Belfast, University Road, Belfast, BT7 1NN, UK

Nick Higginbottom
  Department of Physics and Astronomy, University of Southampton, Southampton, SO17 1BJ, UK

James Matthews
  Department of Physics, Astrophysics, University of Oxford, Denys Wilkinson Building, Keble Road, Oxford, OX1 3RH, UK

Sam Mangham
  Department of Physics and Astronomy, University of Southampton, Southampton, SO17 1BJ, UK

Edward Parkinson
  Department of Electronics and Computer Science, University of Southampton, Southampton, SO17 1BJ, UK

Mandy Hewitt
  School of Mathematics and Physics, Queen's University Belfast, University Road, Belfast, BT7 1NN, UK

Nicolas Scepi
  Univ. Grenoble Alpes, CNRS, IPAG, 38000 Grenoble, France

Austen Wallis
  Department of Physics and Astronomy, University of Southampton, Southampton, SO17 1BJ, UK

Amin Mosallanezhad
  Department of Physics and Astronomy, University of Southampton, Southampton, SO17 1BJ, UK

----------------------------------------

.. toctree::
   :titlesonly:
   :glob:
   :hidden:
   :caption: Documentation

   quick
   installation
   running_sirocco
   input
   output
   plotting
   operation
   radiation
   wind_models
   coordinate
   examples
   physics
   atomic
   meta
   developer
   *
