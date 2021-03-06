% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mfGeometry.R
\name{mfGeomPlanarSizeShape}
\alias{mfGeomPlanarSizeShape}
\title{Planar Size and Shape Space Geometry}
\description{
Geometry of the 'classic' size-and-shape space identifying planar
configurations with a centered and complex vector and treating them as rotation
invariant.
}
\section{Super classes}{
\code{\link[manifoldboost:mfGeometry]{manifoldboost::mfGeometry}} -> \code{\link[manifoldboost:mfGeomEuclidean]{manifoldboost::mfGeomEuclidean}} -> \code{\link[manifoldboost:mfGeomUnitSphere]{manifoldboost::mfGeomUnitSphere}} -> \code{\link[manifoldboost:mfGeomPlanarShape]{manifoldboost::mfGeomPlanarShape}} -> \code{mfGeomPlanarSizeShape}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-mfGeomPlanarSizeShape-register}{\code{mfGeomPlanarSizeShape$register()}}
\item \href{#method-mfGeomPlanarSizeShape-register_v}{\code{mfGeomPlanarSizeShape$register_v()}}
\item \href{#method-mfGeomPlanarSizeShape-log}{\code{mfGeomPlanarSizeShape$log()}}
\item \href{#method-mfGeomPlanarSizeShape-exp}{\code{mfGeomPlanarSizeShape$exp()}}
\item \href{#method-mfGeomPlanarSizeShape-transport}{\code{mfGeomPlanarSizeShape$transport()}}
\item \href{#method-mfGeomPlanarSizeShape-plot}{\code{mfGeomPlanarSizeShape$plot()}}
\item \href{#method-mfGeomPlanarSizeShape-get_normal}{\code{mfGeomPlanarSizeShape$get_normal()}}
\item \href{#method-mfGeomPlanarSizeShape-clone}{\code{mfGeomPlanarSizeShape$clone()}}
}
}
\if{html}{\out{
<details><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="manifoldboost" data-topic="mfGeometry" data-id="distance"><a href='../../manifoldboost/html/mfGeometry.html#method-mfGeometry-distance'><code>manifoldboost::mfGeometry$distance()</code></a></li>
<li><span class="pkg-link" data-pkg="manifoldboost" data-topic="mfGeomEuclidean" data-id="align"><a href='../../manifoldboost/html/mfGeomEuclidean.html#method-mfGeomEuclidean-align'><code>manifoldboost::mfGeomEuclidean$align()</code></a></li>
<li><span class="pkg-link" data-pkg="manifoldboost" data-topic="mfGeomEuclidean" data-id="innerprod"><a href='../../manifoldboost/html/mfGeomEuclidean.html#method-mfGeomEuclidean-innerprod'><code>manifoldboost::mfGeomEuclidean$innerprod()</code></a></li>
<li><span class="pkg-link" data-pkg="manifoldboost" data-topic="mfGeomPlanarShape" data-id="initialize"><a href='../../manifoldboost/html/mfGeomPlanarShape.html#method-mfGeomPlanarShape-initialize'><code>manifoldboost::mfGeomPlanarShape$initialize()</code></a></li>
<li><span class="pkg-link" data-pkg="manifoldboost" data-topic="mfGeomPlanarShape" data-id="structure"><a href='../../manifoldboost/html/mfGeomPlanarShape.html#method-mfGeomPlanarShape-structure'><code>manifoldboost::mfGeomPlanarShape$structure()</code></a></li>
<li><span class="pkg-link" data-pkg="manifoldboost" data-topic="mfGeomPlanarShape" data-id="structure_weights"><a href='../../manifoldboost/html/mfGeomPlanarShape.html#method-mfGeomPlanarShape-structure_weights'><code>manifoldboost::mfGeomPlanarShape$structure_weights()</code></a></li>
<li><span class="pkg-link" data-pkg="manifoldboost" data-topic="mfGeomPlanarShape" data-id="unstructure"><a href='../../manifoldboost/html/mfGeomPlanarShape.html#method-mfGeomPlanarShape-unstructure'><code>manifoldboost::mfGeomPlanarShape$unstructure()</code></a></li>
<li><span class="pkg-link" data-pkg="manifoldboost" data-topic="mfGeomPlanarShape" data-id="unstructure_weights"><a href='../../manifoldboost/html/mfGeomPlanarShape.html#method-mfGeomPlanarShape-unstructure_weights'><code>manifoldboost::mfGeomPlanarShape$unstructure_weights()</code></a></li>
<li><span class="pkg-link" data-pkg="manifoldboost" data-topic="mfGeomPlanarShape" data-id="validate"><a href='../../manifoldboost/html/mfGeomPlanarShape.html#method-mfGeomPlanarShape-validate'><code>manifoldboost::mfGeomPlanarShape$validate()</code></a></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mfGeomPlanarSizeShape-register"></a>}}
\if{latex}{\out{\hypertarget{method-mfGeomPlanarSizeShape-register}{}}}
\subsection{Method \code{register()}}{
center \code{y_}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mfGeomPlanarSizeShape$register(y_)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{y_}}{a complex vector on the sphere, orthogonal to the constant.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mfGeomPlanarSizeShape-register_v"></a>}}
\if{latex}{\out{\hypertarget{method-mfGeomPlanarSizeShape-register_v}{}}}
\subsection{Method \code{register_v()}}{
orthogonally project \code{v_} into the tangent space of \code{y0_}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mfGeomPlanarSizeShape$register_v(v_, y0_ = self$pole_)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{v_}}{a complex tangent vector.}

\item{\code{y0_}}{a complex vector on the sphere, orthogonal to the constant.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mfGeomPlanarSizeShape-log"></a>}}
\if{latex}{\out{\hypertarget{method-mfGeomPlanarSizeShape-log}{}}}
\subsection{Method \code{log()}}{
\code{y_-y0_} after aligning \code{y_} to \code{y0_}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mfGeomPlanarSizeShape$log(y_, y0_ = self$pole_)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{y_}}{a complex vector on the sphere, orthogonal to the constant.}

\item{\code{y0_}}{a complex vector on the sphere, orthogonal to the constant.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mfGeomPlanarSizeShape-exp"></a>}}
\if{latex}{\out{\hypertarget{method-mfGeomPlanarSizeShape-exp}{}}}
\subsection{Method \code{exp()}}{
\code{y0_+v_} assuming \code{v_} in the tangent space of \code{y0_}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mfGeomPlanarSizeShape$exp(v_, y0_ = private$.pole_)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{v_}}{a complex tangent vector.}

\item{\code{y0_}}{a complex vector on the sphere, orthogonal to the constant.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mfGeomPlanarSizeShape-transport"></a>}}
\if{latex}{\out{\hypertarget{method-mfGeomPlanarSizeShape-transport}{}}}
\subsection{Method \code{transport()}}{
size-and-shape parallel transport
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mfGeomPlanarSizeShape$transport(v0_, y0_, y1_, method = c("horizontal"))}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{v0_}}{a complex tangent vector.}

\item{\code{y0_}}{a complex vector on the sphere, orthogonal to the constant.}

\item{\code{y1_}}{a complex vector on the sphere, orthogonal to the constant.}

\item{\code{method}}{currently, only "horizontal" assuming \code{y0_, y1_} aligned
and \code{v0_} a proper horizontal tangent vector at \code{y0_}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mfGeomPlanarSizeShape-plot"></a>}}
\if{latex}{\out{\hypertarget{method-mfGeomPlanarSizeShape-plot}{}}}
\subsection{Method \code{plot()}}{
default plotting function for planar shapes, plotting \code{y_}
in front of \code{y0_} (after alignment).
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mfGeomPlanarSizeShape$plot(
  y_ = self$y_,
  y0_ = self$pole_,
  yaxt = "s",
  bty = "n",
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{y_}}{a complex vector on the sphere, orthogonal to the constant.}

\item{\code{y0_}}{a complex vector on the sphere, orthogonal to the constant.}

\item{\code{...}}{other arguments passed to \code{base::plot}.}

\item{\code{col, pch, type}}{graphical parameters passed to \code{base::plot} referring to \code{y_}.}

\item{\code{ylab, xlab, xlim, ylim, xaxt, yaxt, asp, bty}}{graphical parameters passed to \code{base::plot}
with modified defaults.}

\item{\code{y0_par}}{graphical parameters for \code{y0_}.}

\item{\code{seg_par}}{graphical parameters for line segments connecting \code{y_}
and \code{y0_}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mfGeomPlanarSizeShape-get_normal"></a>}}
\if{latex}{\out{\hypertarget{method-mfGeomPlanarSizeShape-get_normal}{}}}
\subsection{Method \code{get_normal()}}{
Obtain "design matrix" of tangent space normal vectors in
unstructured long format.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mfGeomPlanarSizeShape$get_normal(y0_ = self$pole_, weighted = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{y0_}}{a complex vector on the sphere, orthogonal to the constant.}

\item{\code{weighted}}{logical, should inner product weights be pre-multiplied to
normal vectors?}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mfGeomPlanarSizeShape-clone"></a>}}
\if{latex}{\out{\hypertarget{method-mfGeomPlanarSizeShape-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mfGeomPlanarSizeShape$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
