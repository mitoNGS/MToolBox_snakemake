.. _mtoolbox_variant_annotation:

MToolBox-variant-annotation
===========================

This wrapper performs functional annotation of variants reported in the VCF file (the output of the :ref:`mtoolbox_variant_calling` workflow) with `mtoolnote`_. If the VCF file is not present, the wrapper will first run the :ref:`mtoolbox_variant_calling` workflow to produce it. The final output is an annotated VCF file.

.. note::  If you already have a VCF of mt variants, you might consider to annotate it by directly running `mtoolnote`_.

The setup of this workflow is detailed in :ref:`the setup of the MToolBox-variant-calling workflow<setup_working_directory>`.

.. _`mtoolnote`: https://github.com/mitoNGS/mtoolnote