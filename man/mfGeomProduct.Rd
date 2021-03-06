% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mfGeometry.R
\name{mfGeomProduct}
\alias{mfGeomProduct}
\title{A general product geometry}
\description{
A general product geometry

A general product geometry
}
\details{
An \code{mfGeometry} class implementing a product geometry in an abstract
general way. While there might be more efficient implementations using the
particular structure of your geometry of interest, this should work for
products of any mfGeometries.
}
\section{Super class}{
\code{\link[manifoldboost:mfGeometry]{manifoldboost::mfGeometry}} -> \code{mfGeomProduct}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{mfGeom_default}}{an \code{mfGeometry} object implementing the default
geometry of a component of the product space.}
}
\if{html}{\out{</div>}}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{y_}}{return lists of respective elements in the components.
\code{y_} can be specified as a list of \code{mfGeometry} objects constituting
the components of the product geometry.}

\item{\code{pole_}}{return lists of respective elements in the components.}

\item{\code{weights_}}{return lists of respective elements in the components.}

\item{\code{mfGeom_classes}}{return a list of the classes of \code{y_} (read only).}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-mfGeomProduct-new}{\code{mfGeomProduct$new()}}
\item \href{#method-mfGeomProduct-slice}{\code{mfGeomProduct$slice()}}
\item \href{#method-mfGeomProduct-structure}{\code{mfGeomProduct$structure()}}
\item \href{#method-mfGeomProduct-unstructure}{\code{mfGeomProduct$unstructure()}}
\item \href{#method-mfGeomProduct-unstructure_weights}{\code{mfGeomProduct$unstructure_weights()}}
\item \href{#method-mfGeomProduct-structure_weights}{\code{mfGeomProduct$structure_weights()}}
\item \href{#method-mfGeomProduct-align}{\code{mfGeomProduct$align()}}
\item \href{#method-mfGeomProduct-register}{\code{mfGeomProduct$register()}}
\item \href{#method-mfGeomProduct-register_v}{\code{mfGeomProduct$register_v()}}
\item \href{#method-mfGeomProduct-distance}{\code{mfGeomProduct$distance()}}
\item \href{#method-mfGeomProduct-exp}{\code{mfGeomProduct$exp()}}
\item \href{#method-mfGeomProduct-log}{\code{mfGeomProduct$log()}}
\item \href{#method-mfGeomProduct-transport}{\code{mfGeomProduct$transport()}}
\item \href{#method-mfGeomProduct-innerprod}{\code{mfGeomProduct$innerprod()}}
\item \href{#method-mfGeomProduct-plot}{\code{mfGeomProduct$plot()}}
\item \href{#method-mfGeomProduct-get_normal}{\code{mfGeomProduct$get_normal()}}
\item \href{#method-mfGeomProduct-validate}{\code{mfGeomProduct$validate()}}
\item \href{#method-mfGeomProduct-clone}{\code{mfGeomProduct$clone()}}
}
}
\if{html}{\out{
<details open><summary>Inherited methods</summary>
<ul>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mfGeomProduct-new"></a>}}
\if{latex}{\out{\hypertarget{method-mfGeomProduct-new}{}}}
\subsection{Method \code{new()}}{
Initialize product geometry with each
component the \code{mfGeom_default} based on the
\code{data} containing the variables for \code{y_}
specified in the \code{formula}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mfGeomProduct$new(mfGeom_default, data, formula)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{mfGeom_default}}{\code{mfGeometry} object carrying the
component geometry.}

\item{\code{data}}{containing required variables.}

\item{\code{formula}}{formula containing interpreted
via \code{mfInterpret_formula}) and indicating, in particular,
the ID variable for splitting observations in the product
components.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mfGeomProduct-slice"></a>}}
\if{latex}{\out{\hypertarget{method-mfGeomProduct-slice}{}}}
\subsection{Method \code{slice()}}{
subset product geometry to specified components.
Make sure to clone the product geometry first if
complete geometry should be preserved.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mfGeomProduct$slice(which)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{which}}{integers or character strings indicating
which elements of \code{y_} should be selected.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mfGeomProduct-structure"></a>}}
\if{latex}{\out{\hypertarget{method-mfGeomProduct-structure}{}}}
\subsection{Method \code{structure()}}{
bring \code{y} into the right order to
split and pass it to the component geometries.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mfGeomProduct$structure(y)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{y}}{a list of \code{y} elements of the component manifolds in an
unstructured 'flattened' way}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mfGeomProduct-unstructure"></a>}}
\if{latex}{\out{\hypertarget{method-mfGeomProduct-unstructure}{}}}
\subsection{Method \code{unstructure()}}{
unstructure \code{y_} in the component geometries
and arrange it into the right order.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mfGeomProduct$unstructure(y_)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{y_}}{a list of objects living in the component manifolds}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mfGeomProduct-unstructure_weights"></a>}}
\if{latex}{\out{\hypertarget{method-mfGeomProduct-unstructure_weights}{}}}
\subsection{Method \code{unstructure_weights()}}{
unstructure \code{weights_} in the component geometries
and arrange it into the right order.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mfGeomProduct$unstructure_weights(weights_)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{weights_}}{list of numeric vectors of inner product weights.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mfGeomProduct-structure_weights"></a>}}
\if{latex}{\out{\hypertarget{method-mfGeomProduct-structure_weights}{}}}
\subsection{Method \code{structure_weights()}}{
bring \code{weights} into the right order to
split and pass it to the component geometries.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mfGeomProduct$structure_weights(weights)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{weights}}{numeric vector of inner product weights matching the \code{y}
as a numeric vector.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mfGeomProduct-align"></a>}}
\if{latex}{\out{\hypertarget{method-mfGeomProduct-align}{}}}
\subsection{Method \code{align()}}{
component-wise alignment.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mfGeomProduct$align(y_, y0_, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{y_}}{a list of objects living in the component manifolds}

\item{\code{y0_}}{a list of objects living in the component manifolds}

\item{\code{...}}{other arguments passed to underlying \code{mfGeometry}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mfGeomProduct-register"></a>}}
\if{latex}{\out{\hypertarget{method-mfGeomProduct-register}{}}}
\subsection{Method \code{register()}}{
component-wise registration.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mfGeomProduct$register(y_, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{y_}}{a list of objects living in the component manifolds}

\item{\code{...}}{other arguments passed to underlying \code{mfGeometry}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mfGeomProduct-register_v"></a>}}
\if{latex}{\out{\hypertarget{method-mfGeomProduct-register_v}{}}}
\subsection{Method \code{register_v()}}{
component-wise registration.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mfGeomProduct$register_v(v_, y0_, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{v_}}{a list of tangent vectors living in the component manifolds}

\item{\code{y0_}}{a list of objects living in the component manifolds}

\item{\code{...}}{other arguments passed to underlying \code{mfGeometry}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mfGeomProduct-distance"></a>}}
\if{latex}{\out{\hypertarget{method-mfGeomProduct-distance}{}}}
\subsection{Method \code{distance()}}{
compute vector of distances in component
geometries.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mfGeomProduct$distance(y0_, y1_, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{y0_}}{a list of objects living in the component manifolds}

\item{\code{y1_}}{a list of objects living in the component manifolds}

\item{\code{...}}{other arguments passed to underlying \code{mfGeometry}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mfGeomProduct-exp"></a>}}
\if{latex}{\out{\hypertarget{method-mfGeomProduct-exp}{}}}
\subsection{Method \code{exp()}}{
loop over Exp-maps in the component
geometries.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mfGeomProduct$exp(v_, y0_, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{v_}}{a list of tangent vectors living in the component manifolds}

\item{\code{y0_}}{a list of objects living in the component manifolds}

\item{\code{...}}{other arguments passed to underlying \code{mfGeometry}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mfGeomProduct-log"></a>}}
\if{latex}{\out{\hypertarget{method-mfGeomProduct-log}{}}}
\subsection{Method \code{log()}}{
loop over Log-maps in the component
geometries.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mfGeomProduct$log(y_, y0_, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{y_}}{a list of objects living in the component manifolds}

\item{\code{y0_}}{a list of objects living in the component manifolds}

\item{\code{...}}{other arguments passed to underlying \code{mfGeometry}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mfGeomProduct-transport"></a>}}
\if{latex}{\out{\hypertarget{method-mfGeomProduct-transport}{}}}
\subsection{Method \code{transport()}}{
loop over parallel transports in the component
geometries.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mfGeomProduct$transport(v0_, y0_, y1_, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{v0_}}{a list of tangent vectors living in the component manifolds}

\item{\code{y0_}}{a list of objects living in the component manifolds}

\item{\code{y1_}}{a list of objects living in the component manifolds}

\item{\code{...}}{other arguments passed to underlying \code{mfGeometry}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mfGeomProduct-innerprod"></a>}}
\if{latex}{\out{\hypertarget{method-mfGeomProduct-innerprod}{}}}
\subsection{Method \code{innerprod()}}{
compute vector of inner products in the
component geometries. An inner product on the product
space can be obtained as a scalar product of the
returned vector.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mfGeomProduct$innerprod(v0_, v1_ = v0_, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{v0_}}{a list of tangent vectors living in the component manifolds}

\item{\code{v1_}}{a list of tangent vectors living in the component manifolds}

\item{\code{...}}{other arguments passed to underlying \code{mfGeometry}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mfGeomProduct-plot"></a>}}
\if{latex}{\out{\hypertarget{method-mfGeomProduct-plot}{}}}
\subsection{Method \code{plot()}}{
loop over individual plots of the components.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mfGeomProduct$plot(y_ = self$y_, y0_ = self$pole_, main = names(y_), ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{y_}}{a list of objects living in the component manifolds}

\item{\code{y0_}}{a list of objects living in the component manifolds}

\item{\code{main}}{vector of plot titles}

\item{\code{...}}{other arguments passed to underlying \code{mfGeometry}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mfGeomProduct-get_normal"></a>}}
\if{latex}{\out{\hypertarget{method-mfGeomProduct-get_normal}{}}}
\subsection{Method \code{get_normal()}}{
Obtain "design matrix" of tangent space normal vectors in
unstructured long format, arranging rows into the right
order after evaluating it on the components.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mfGeomProduct$get_normal(y0_ = self$pole_, weighted = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{y0_}}{a list of objects living in the component manifolds}

\item{\code{weighted}}{logical, should inner product weights be pre-multiplied to
normal vectors?}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mfGeomProduct-validate"></a>}}
\if{latex}{\out{\hypertarget{method-mfGeomProduct-validate}{}}}
\subsection{Method \code{validate()}}{
loop of validation functions of component geometries.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mfGeomProduct$validate(y_)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{y_}}{a list of objects living in the component manifolds}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mfGeomProduct-clone"></a>}}
\if{latex}{\out{\hypertarget{method-mfGeomProduct-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mfGeomProduct$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
