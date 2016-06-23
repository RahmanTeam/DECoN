\name{NEWS}
\title{News for \R Package \pkg{Matrix}}% MM: look into ../svn-log-from.all
\encoding{UTF-8}
\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}
%% NB: The date (yyyy-mm-dd) is the "Packaged:" date in ../DESCRIPTION

\section{Changes in version 1.2-0 (2015-04-03)}{
  \subsection{New Features}{
    \itemize{
      \item New \code{\%&\%} for \dQuote{boolean arithmetic} matrix product.
      \item New argument \code{boolArith = NA} in \code{crossprod()} and
      \code{tcrossprod()}.  \code{boolArith = TRUE} now forces boolean
      arithmetic, where \code{boolArith = FALSE} forces numeric one.

      Several of these products are more efficient thanks to new C
      functionality based on our new \code{chm_transpose_dense()}, and
      others based on \code{geMatrix_crossprod},
      \code{geMatrix_matrix_mm}, etc.

      \item Most dense matrix products, also for non-\code{dgeMatrix},
      including \code{"l..Matrix"} and \code{"n..Matrix"} ones are now
      directly handled by new \code{.Call()}s.

      \item \code{"dMatrix"} (numeric) and \code{"lMatrix"} (logical)
      matrices can now be coerced to \code{"nMatrix"} (non-zero pattern
      or \dQuote{boolean}) even when they contain \code{NA}s, which then
      become \code{TRUE}s.

      \item More thorough checking of \code{cbind2()} and
      \code{rbind2()} methods, notably as they are called from \code{cbind()}
      and \code{rbind()} from \R version 3.2.0 on.

      \code{rbind2(<dense>, <dense>)} is faster, being based on new C code.

      \item symmetric Matrices (i.e., inheriting from
      \code{"symmetricMatrix"}) are allowed to have \code{dimnames} of
      the form \code{list(NULL, <names>)} \emph{and} now print correctly
      and get correctly coerced to general matrices.

      \item \code{indMatrix} object (\dQuote{index matrices}) no longer
      need to be \dQuote{skinny}.

      \item \code{rsparseMatrix()} now accepts \code{rand.x = NULL} and
      then creates a random \emph{patter\bold{n}} matrix
      (\code{"nsparseMatrix"}).

      \item \code{anyDuplicatedT()} and \code{uniqTsparse()} low level
      utilities are exported now.

      \item Partial Korean translations of messages.
    }
  }
  \subsection{Deprecation}{
    \itemize{
      \item For \eqn{R \ge 3.2.0}, \code{cBind()} and \code{rBind()}
      are deprecated, as they are no longer needed since \code{cbind()}
      and \code{rbind()} do work automatically.
    }
  }

  \subsection{Bug Fixes}{
    \itemize{
      \item Fix some \code{rbind2()} methods.
      \item \code{t()} now tranposes the dimnames even for symmetric
      matrices.
      \item \code{diag(M) <- val} did not always recycle \code{val} to
      full length, e.g., when \code{M} was a \code{"dtrMatrix"}.
      \item \code{crossprod(<indMatrix>)} was wrong in cases where the
      matrix had all-zero columns.
      \item Matrix products (\code{\%*\%}, \code{crossprod}, and
      \code{tcrossprod}) with one sparse and one dense argument now
      return \emph{numeric} (a \code{"dMatrix"}) when they should, i.e.,
      unless the new setting \code{boolArith = TRUE} is applied.
    }
  }
}

\section{Changes in version 1.1-5 (2015-01-18, svn r3037)}{
  \subsection{New Features}{
    \itemize{
      \item Start using \code{anyNA()} (for speedup).
      \item Matrix products (\code{\%*\%}, \code{crossprod},
      \code{tcrossprod}) now behave compatibly to \R 3.2.0, i.e., more
      lenient in matching dimensions for matrix - vector products.
      \item \code{isTriangular()} gets new optional argument \code{upper = NA}.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item \code{crossprod()} and \code{tcrossprod()} fixes for
      several <diagonal> o <sparse> combinations.
      \item \code{rowMeans(<dgeMatrix>, na.rm=TRUE)} was wrong sometimes.
      \item fix and speedup of coercions (\code{as(., .)}) from and to
      symmetric or triangular matrices.
      \item \code{invPerm()} coercion to integer
      \item \code{dimnames( solve(.,.) )} fix [r3036]
      \item \code{tril()} and \code{triu()} now return correct \code{uplo}.
      \item \code{names(dimnames(.))} now preserved, e.g. in
      \code{symmpart()} or subsetting (\code{A[i,j]}).
    }
  }
}

\section{Changes in version 1.1-4 (2014-06-14, svn r2994)}{
  \subsection{New Features}{
    \itemize{
      \item new \code{rsparsematrix()} for random sparse Matrices.
      \item improved warnings, notably for unused arguments previously
      swallowed into \code{...}.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item \code{crossprod(<vec>, <dsyMatrix>)} fixed.
      \item \code{crossprod()} and \code{kronecker()} fixes for some
      <indMatrix> cases.
    }
  }
}

\section{Changes in version 1.1-3 (2014-03-30, svn r2982)}{
  \subsection{New Features}{
    \itemize{
      \item \code{\%*\%} and \code{crossprod()} now also work with
      \code{sparseVector}s.
      \item speedup of \code{crossprod(v, <sparseM>)}, thanks to nudge
      by Niels Richard Hansen.
      \item new help page for all such matrix products
      (\file{../man/matrix-products.Rd}).
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item \code{image()} now get correct \code{ylim} again.
      \item More consistent matrix products.
    }
  }
}

\section{Changes in version 1.1-2-2 (2014-03-04, svn r2966)}{
  \subsection{Bug Fixes}{
    \itemize{
      \item correct adaption to \R 3.1.0
      \item using \code{tolerance} (and not \sQuote{tol}) in \code{all.equal()}
    }
  }
}

\section{Changes in version 1.1-2 (2014-01-28, svn r2962)}{
  \subsection{New Features}{
    \itemize{
      \item export fast power-user coercion utilities
      \code{.dsy2mat()}, \code{.dxC2mat()}, \code{.T2Cmat()}, \code{..2dge()}.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item matrix products now (mostly) work with \code{sparseVector}s;
      and correctly in some more cases.
    }
  }
}

\section{Changes in version 1.1-1.1 (2013-12-30, svn r2957)}{
    \itemize{
      \item Testing code's \code{assertWarning()} adapted for \eqn{R \le 3.0.1}.
      \item \code{Depends: R >= 2.15.2} eases checking.
    }
}

\section{Changes in version 1.1-1 (2013-12-28)}{
  \subsection{New Features}{
    \itemize{
      \item \code{image(.., xlim, ylim)}: nicer defaults %% ../R/dgTMatrix.R
      for the axis limits, and \code{ylim} is sorted decreasingly; not
      strictly back-compatible but should never harm.
      \item \code{rankMatrix(*, method="qr")} now using \code{tol}
      \item \code{T2graph()} and \code{graph2T()} export old functionality explicitly.
      Tweaks in conversions between \code{"graph"} and
      \code{"sparseMatrix"} objects.  Notably, \code{as(<graph>,
	<Matrix>)} now more often returns a (0/1 pattern) "n..Matrix".
      \item \code{sparseMatrix()}: new \code{use.last.ij} argument.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item \code{KhatriRao()}: fix rownames (X <-> Y)
      \item \code{qr.coef()}, \code{qr.fitted}, and \code{qr.resid} now
      also work with \emph{sparse} RHS \code{y}.
      \item sparse matrix \dQuote{sub assignments}, e.g., \code{M[ii] <- v},
      speedup and fixes.
      \item bug fixes also in \code{M[negative indices] <- value} and
      \code{<sparseMatrix>[cbind(i,j)]}.
    }
  }
}

\section{Changes in version 1.1-0 (2013-10-21, svn r2930)}{
  \subsection{New Features}{
    \itemize{
      \item \code{fac2sparse} and \code{fac2Sparse} now exported, with a
      new \code{giveCsparse} option.
      \item Update to latest \command{SuiteSparse} C library by Tim Davis,
      U. Florida.
      \item ensuing \dQuote{C API changes}
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item fixed long lasting undetected \code{solve(<dsCMatrix>, *)} bug
      \item
    }
  }
}


\section{Changes in version 1.0-14 (2013-09-12)}{
  \subsection{New Features}{
    \itemize{
      \item .
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item .
    }
  }
}

\section{Changes in version 1.0-13 (2013-09-10)}{
  \subsection{New Features}{
    \itemize{
      \item New efficient \code{KhatriRao()} by Michael Cysouw
      \item
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item .
      \item
    }
  }
}

\section{Changes in version 1.0-12 (2013-03-26)}{
  \subsection{New Features}{
    \itemize{
      \item .
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item .
      \item
    }
  }
}


\section{Changes in version 1.0-11 (2013-02-02)}{
  \subsection{New Features}{
    \itemize{
      \item .
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item \code{as(<csr>, "dgCMatrix")} (from package \pkg{SparseM})
      now works again.
      \item .
    }
  }
}

\section{Changes in version 1.0-10 (2012-10-22)}{
  \subsection{New Features}{
    \itemize{
      \item \code{.sparseDiagonal()}: new \code{unitri} argument, and
      more flexibility;
      \item new \code{solve(<dsCMatrix>, <missing>)} via efficient C code.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item .
    }
  }
}

\section{Changes in version 1.0-9 (2012-09-05)}{
  \subsection{New Features}{
    \itemize{
      \item new \code{sparseVector()} constructor function.
      \item \code{is.finite()} \code{is.infinite()} now work for our
      matrices and "*sparseVector" objects.
      \item \code{diag(.) <- V} now preserves symmetricity,
      triangularity and even uni-triangularity sometimes.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item Quite a few fixes for \code{Ops} (arithmetic, logic, etc)
      group methods.
      \item Ditto for \code{diagonalMatrix} methods.
    }
  }
}

\section{Changes in version 1.0-6 (2012-03-16, publ. 2012-06-18)}{
  \subsection{New Features}{
    \itemize{
      \item .
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item .
    }
  }
}

\section{Changes in version 1.0-5 (2012-03-15)}{
  \subsection{New Features}{
    \itemize{
      \item .
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item .
    }
  }
}

\section{Changes in version 1.0-4 (2012-02-21)}{
  \subsection{New Features}{
    \itemize{
      \item .
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item .
    }
  }
}

\section{Changes in version 1.0-3 (2012-01-13)}{
  \subsection{New Features}{
    \itemize{
      \item .
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item .
    }
  }
}

\section{Changes in version 1.0-2 (2011-11-19)}{
  \subsection{New Features}{
    \itemize{
      \item .
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item .
    }
  }
}

\section{Changes in version 1.0-1 (2011-10-18)}{
  \subsection{New Features}{
    \itemize{
      \item .
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item .
    }
  }
}

\section{Changes in version 1.0-0 (2011-10-04)}{
  \subsection{New Features}{
    \itemize{
      \item .
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item .
    }
  }
}

\section{Changes in version 0.9996875-3 (2011-08-13)}{
  \subsection{New Features}{
    \itemize{
      \item .
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item .
    }
  }
}

\section{Changes in version 0.9996875-2 (2011-08-09)}{
  \subsection{New Features}{
    \itemize{
      \item .
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item .
    }
  }
}

\section{Changes in version 0.9996875-1 (2011-08-08)}{
  \subsection{New Features}{
    \itemize{
      \item .
    }
  }

  \subsection{Bug Fixes}{
    \itemize{
      \item .
    }
  }
}

\section{Changes in version 0.999375-50 (2011-04-08)}{
  \subsection{New Features}{
    \itemize{
      \item .
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item .
    }
  }
}
% How can I add vertical space ?
% \preformatted{} is not allowed, nor is \cr

%--------------- pre-pre-history: ------------------------


\section{Version 0.2-1 (2000-07-15)}{
  The first CRAN release of the \pkg{Matrix} package,
  titled \dQuote{A Matrix library for R} authored by
  Douglas Bates (maintainer, principal author) and Saikat DebRoy.
  \subsection{Features}{
    \itemize{
      \item \code{Matrix()} constructor for \R objects of class \code{Matrix}.
      \item \code{Matrix.class()} returning informal subclasses such as
      \code{"Hermitian"}, \code{"LowerTriangular"}
      \item \code{is.Orthonormal()}, \code{is.Hermitian()} ,
      \code{is.UpperTriangular()} functions.
      \item \code{SVD()}, \code{lu()}, and \code{schur()} decomposition
      generics with \code{"Matrix"} methods.
      \item \code{rcond()}, \code{norm()}, \code{det()};
      \code{expand()} and \code{facmul()}.

      \item C++ interface to LAPACK
    }
  }
}
