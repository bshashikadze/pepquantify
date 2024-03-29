% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_diann.R
\name{read_diann}
\alias{read_diann}
\title{Read and filter DIA precursor-level data produced by the DIA-NN}
\usage{
read_diann(
  Q_Val = 0.01,
  Global_Q_Val = 0.01,
  Global_PG_Q_Val = 0.01,
  Lib_Q_Val = 0.01,
  Lib_PG_Q_Val = 0.01,
  experimental_library,
  unique_peptides_only = TRUE,
  Quant_Qual = 0.5,
  diann_file_name = "",
  directory = ".",
  id_column = "Genes",
  second_id_column = "Protein.Group",
  quantity_column = "Genes.MaxLFQ.Unique",
  for_msempire = T,
  sum_charge = TRUE,
  save_supplementary = TRUE,
  exclude_samples = c(),
  include_mod_in_pepreport = F
)
}
\arguments{
\item{Q_Val}{refer to https://github.com/vdemichev/DiaNN}

\item{Global_Q_Val}{refer to https://github.com/vdemichev/DiaNN}

\item{Global_PG_Q_Val}{refer to https://github.com/vdemichev/DiaNN}

\item{Lib_Q_Val}{refer to https://github.com/vdemichev/DiaNN}

\item{Lib_PG_Q_Val}{refer to https://github.com/vdemichev/DiaNN}

\item{experimental_library}{set true if you use empirical libraries (e.g. prefractionation or GPF), false in case of lib free search with mbr enabled}

\item{unique_peptides_only}{TRUE only unique peptides will be used for quantification (recommended)}

\item{Quant_Qual}{refer to https://github.com/vdemichev/DiaNN; pepquantify by default sets it to 0.5}

\item{diann_file_name}{by default function will assume that largest .tsv file is the DIA-output, however specific file name can be specified diann_file_name=filename.tsv}

\item{directory}{default is a current directory}

\item{id_column}{default "Genes" (can be switched to Protein.Group)}

\item{second_id_column}{default "Protein.Group" (can be switched to Genes)}

\item{quantity_column}{default "Genes.MaxLFQ.Unique", not important for MS-EmpiRe}

\item{for_msempire}{default TRUE so pepquantify will prepare data for MS-EmpiRe quantification, if purpose is to filter DIA-NN output and generate peptides and protein groups file set to false}

\item{sum_charge}{how precursor charge states will be aggregated to peptide level, True means the sum will be taken, in case of false, precursor with the highest intensity will be kept, default false}

\item{save_supplementary}{default TRUE, output is peptide and protein level data which can be used as a supplement}

\item{exclude_samples}{if not empty, excludes specified sample/s from further analysis (only if necessary, e.g. after inspecting PCA)}

\item{include_mod_in_pepreport}{default FALSE, if true includes modifications in the output peptide file (currently only Carbamidomethyl (C))}
}
\value{
The list of three elements, the first is the filtered peptides file, the second is the filtered protein groups file and the last is the character that stores the type of experiment
}
\description{
Read and filter DIA precursor-level data produced by the DIA-NN
}
\examples{
read_diann(exclude_samples=c("samplename"), experimental_library = TRUE)
}
