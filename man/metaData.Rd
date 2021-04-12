\name{metaData}
\docType{data}
\alias{metaData}
\title{Dose response data from several published meta-analyses}
\description{
Dose response data from over 200 compounds included in published meta-analyses.  The data are aggregated in a single data frame in a common format.
}
\usage{data('metaData')
}
\format{
The data frame has one row for each compound, protocol within compound, and 
dose group within protocol.    
Compound and protocol level descriptors are repeated on each
row of the data frame.
  \describe{
		\item{\code{drugid}}{A numerical ID identifying each drug}
		\item{\code{taid}}{A drug can be studied in more than one therapeutic area.  
		The \code{taid} ID identifies each TA/drug combination.} 
		\item{\code{protid}}{Numerical (1,2,3,...) ID for protocols specific to each TAID.}
		\item{\code{gname}}{Generic drug name}
		\item{\code{bname}}{Branded(USA) drug name}
		\item{\code{drugtype}}{Drug classified as SMALL MOLECULE, BIOLOGIC, OTHER}
		\item{\code{route}}{Route of administration, e.g., oral, subcutaneous,...}
		\item{\code{routeShort}}{Abbreviated format for \code{route}}
		\item{\code{oralForm}}{Formulation (e.g., TABLET, POWDER,...) for drugs with oral administration.}
		\item{\code{fdaapproved}}{NA if status was not yet determined}
		\item{\code{metasource}}{Meta-analysis contributing compounds.  BIO14: biological compounds through 2014;
		FDA914: FDA approved small molecules and 'other' 2009-2014;  FDA1417: FDA approved compounds 2014-2017;
		Pfizer P2 compounds 1998-2009; PFIZERUPDATE18:  Pfizer compounds 2009-2018}
		\item{\code{protno}}{Sponsor assigned protocol name/number}
		\item{\code{nctno}}{Clintrial.gov protocol ID}
		\item{\code{protyear}}{When available, year of first patient/first visit.  In some cases, date of journal publication}
		\item{\code{design}}{PARELLEL, CROSSOVER,...}
		\item{\code{actcomp}}{Indicator if an active comparator was included in the protocol}
		\item{\code{etype}}{\code{etype=1} for the designated primary endpoint.  For completeness, where
		there was ambiguity in the selection of the endpoint, additional endpoint data was included on separate 
		rows and indicated by etype=2,3,...  Most analyses subset on \code{etype=1}}
		\item{\code{poptype}}{For a compound and TA, there can be distinctly different populations with 
		anticipated response differences, e.g., treatment-naive and pre-treated patients.  The population with 
		the most studied doses has \code{poptype=1}.  For completeness, additional populations are included 
		and identified by \code{poptype=2,3,...}.  Most analyses subset on \code{poptype=1}}
		\item{\code{primsource}}{IRO/PRO investigator/patient reported outcome; L lab, V vitals}
		\item{\code{primtype}}{Primary endpoint is BINARY, CONTINOUS, TIMETOEVEN}
		\item{\code{primtime}}{time units to primary endpoint from randomization}
		\item{\code{timeunit}}{DAY, HR, MIN, MONTH, WK for primary endpoint}
		\item{\code{indication}}{Disease description}
		\item{\code{broadta}}{Broad TA classification of the indication}
		\item{\code{endpointLong,endpointShort}}{Endpoint name and an abbreviated form using 
		for example, cfb and pcfb for 
		change from baseline and percent change from baseline}
		\item{\code{dose}}{Total daily dose for small molecules, total weekly dose for biologics in mg or mg/kg 
		for weight-based dosing.}
		\item{\code{tload}}{Amount of any loading dose}
		\item{\code{nload}}{Number of visits with a loading dose}
		\item{\code{regimen}}{Dosing frequency}
		\item{\code{primregimen}}{\code{primregimen=1} for most doses/regimens, but \code{primregimen=2} for a few
		regimens that clearly differed from the most common regimen for the same total dose. Most analyses 
		subset on \code{primregimen=1}}
		\item{\code{rslt}}{The sample dose group mean (continuous) or proportion (binary) of the primary endpoint.
		Analyses of the
		time-to-event endpoints was compound specific (either a mean or a proportion was estimated). 
		}
		\item{\code{se}}{Standard error of \code{rslt}}
		\item{\code{sd}}{Dose group sample standard deviation for continuous data}
		\item{\code{lcl, ucl, alpha}}{alpha-level interval (lcl,ucl) when confidence intervals 
		were extracted from the original data source because se were not reported}
		\item{\code{sampsize}}{Sample size reported for \code{rslt}.  The handling of missing data by
		the protocol sponsors varied, but 'completers' was most common.}
		\item{\code{ittsize}}{The number randomized. The counts are usually available, except for 
		internal data before 2009, where it was not collected.}
		\item{\code{pmiss}}{Percent of missing data.}
  }
}

\details{
Compound sampling plans and other details are given in the publications:

Thomas, N., Sweeney, K., and Somayaji, V. (2014).  Meta-analysis of clinical dose response in a large drug development portfolio, Statistics in Biopharmaceutical Research, Vol. 6, No.4, 302-317.
<doi:10.1080/19466315.2014.924876>

Thomas, N., and Roy, D. (2016). Analysis of clinical dose-response in small-molecule drug development: 2009-2014. Statistics in Biopharmaceutical Research, Vol. 6, No.4, 302-317 <doi:10.1080/19466315.2016.1256229>

Wu, J., Banerjee, A., Jin, B. Menon, M. S., Martin, S. and Heatherington, A. (2017). Clinical dose response for a broad set of biological products: A model-based meta-analysis. Statistical Methods in Medical Research. <doi:10.1177/0962280216684528>
}

\examples{
data('metaData')
names(metaData)
}

\keyword{datasets}
