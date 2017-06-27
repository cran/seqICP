%% Define Macros to copy information from DESCRIPTION file
% To get the package title at build time from the DESCRIPTION file
\newcommand{\packageTitle}{\Sexpr[results=rd,stage=build]{tools:::Rd_package_title("#1")}}
% To get the package description at build time from the DESCRIPTION file
\newcommand{\packageDescription}{\Sexpr[results=rd,stage=build]{tools:::Rd_package_description("#1")}}
% To get the package author at build time from the DESCRIPTION file
\newcommand{\packageAuthor}{\Sexpr[results=rd,stage=build]{tools:::Rd_package_author("#1")}}
% To get the package maintainer at build time from the DESCRIPTION file
\newcommand{\packageMaintainer}{\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer("#1")}}
% To get a formatted copy of the whole DESCRIPTION file
\newcommand{\packageDESCRIPTION}{\Sexpr[results=rd,stage=build]{tools:::Rd_package_DESCRIPTION("#1")}}
% To include various indices about an installed package
\newcommand{\packageIndices}{\Sexpr[results=rd,stage=build]{tools:::Rd_package_indices("#1")}}

\name{seqICP-package}

\alias{seqICP-package}

\alias{seqICP_package}

\docType{package}

\title{\packageTitle{seqICP}}

\description{\packageDescription{seqICP}}

\details{
The DESCRIPTION file:
\packageDESCRIPTION{seqICP}
\packageIndices{seqICP}
}

\author{
\packageAuthor{seqICP}

Maintainer: \packageMaintainer{seqICP}
}

\references{
Pfister, N., P. Bühlmann and J. Peters (2017).
Invariant Causal Prediction for Sequential Data. ArXiv e-prints (1706.08058).

Peters, J., P. Bühlmann, and N. Meinshausen (2016).
Causal inference using invariant prediction: identification and confidence intervals.
Journal of the Royal Statistical Society, Series B (with discussion) 78 (5), 947–1012.
}

\keyword{ package }
