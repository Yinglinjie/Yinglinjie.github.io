layout: post
title: "FEM"
date: 2023-02-13 
categories: CATEGORY-1 CATEGORY-2

::: flushleft
Student Name: Linjie Ying\
Instructor: Dr. Pengtao Sun
:::

::: flushright
![image](png/UNLV.png){height="1cm"}
:::

::: center
**Notes on Finite Element Method**\
Assignment 1: Galerkin method on simple elliptic equations
:::

------------------------------------------------------------------------

We intend to build a weak solution of the elliptic problem
$$\label{eq:ellip}
  \begin{cases}
    -\Delta u+au&= f\textmd{ in }\Omega\\
    \alpha u+\beta \frac{\partial u}{\partial n}&=0\textmd{ on }\partial \Omega
  \end{cases}$$ by first constructing solutions of certain
finite-dimensional approximations to a weak form of
([\[eq:ellip\]](#eq:ellip){reference-type="ref" reference="eq:ellip"})
and then passing to limits.

# Weak solution {#sec:weak}

##  Poisson equation with Dirichlet boundary conditions {#sec:dirichlet}

First model problem: $$\label{eq:homo}
  \begin{cases}
    -\Delta u&= f\textmd{ in }\Omega\\
    u&=0\textmd{ on }\partial \Omega
  \end{cases}$$ where $f\in L^2(\Omega)$ and $\Omega$ is open bounded
subset of $\mathbb{R}^n$ with $C^1$ boundary $\partial \Omega$. In order
to define a weak solution, we assume first that $u\in C^\infty_0$, then
multiply the equation $-\Delta u= f$ by a smooth test function
$v\in C_0^\infty$, and integrate over $\Omega$ to find
$$\int_\Omega \nabla u\cdot \nabla v=\int _\Omega fv,$$ where we have
integrated by parts on left side. By approximation
([7](#thm:approx){reference-type="ref" reference="thm:approx"}) we find
the same identity holds with the smooth function $v$ replaced by any
$v\in H^1_0(\Omega)$, and the resulting identity
([\[eq:homo\]](#eq:homo){reference-type="ref" reference="eq:homo"})
makes sense if only $u\in H^1_0(\Omega)$. Then we get a bilinear form
associated with the elliptic operator $-\Delta u$ as
$$\label{eq:variation}
  a(u,v)=\int_\Omega \nabla u\cdot \nabla v$$ for
$u,v\in H^1_0(\Omega)$.

::: defn
**Definition 1**. *We say that $u\in H_0^1(\Omega)$ is a weak solution
of the boundary-value problem $(\ref{eq:homo})$ if $$\label{eq:equal}
    a(u,v)=f(v)\quad \textmd{ for all } v\in H^1_0(\Omega),$$*

*where $(\ ,\ )$ denotes the inner product in $L^2(\Omega),$ and $f(v)$
is a linear form of $f:H_0^1\rightarrow \mathbb{R}$.*
:::

#### Existence and uniqueness of weak solution

Denote semi-norm of $u\in H^m(\Omega)$ as
$$|u|_{m,\Omega}=\left(\sum_{|\alpha|=m} \int_{\Omega}|\partial^\alpha u|^2\right)^{\frac{1}{2}}.$$
We see that $$a(u,u)=|\nabla u|_{0,\Omega}^2$$ Futhermore, $a(u,u)=0$
implies $Du = 0$ and thus $u$ is constant. As $u=0$ on $\partial \Omega$
we have $u\equiv 0$ on $\Omega.$ Hence, $a(\cdot,\cdot )$ defines an
inner product on $H^1_0$, and then the problem
([\[eq:equal\]](#eq:equal){reference-type="ref" reference="eq:equal"})
has a unique solution by Riesz Representation theorem.

Generally, use the Lax-Milgram Theorem, which is primarily significant
when it does not require symmetry of $a(\cdot, \cdot)$. Check the two
condition of $a(u,v)$, it's obvious that
$$|a(u,v)|\leq |\nabla u|_{0,\Omega}|\nabla v|_{0,\Omega}\leq \Vert u\Vert_{H^1_0(\Omega)}\Vert v\Vert_{H^1_0(\Omega)}.$$
Use Poincare's inequality that
$$\Vert u\Vert_{L^2(\Omega)} \leq C \Vert \nabla u\Vert_{L^2(\Omega)},$$
it easily follows that
$$\frac{1}{C}\Vert u\Vert_{H^1_0(\Omega)^2}\leq \Vert \nabla u\Vert_{L^2(\Omega)} = a(u,u).$$
Futhermore the operator $f(v), v\in H_0^1(\Omega)$ is continuous by the
Cauchy-Schwarz inequality, $$\begin{aligned}
    |f(v)|\leq \int_\Omega|f||v|\leq \Vert f\Vert_{L^2(\Omega)}\Vert v\Vert_{L^2(\Omega)}
    \leq \Vert f\Vert_{L^2(\Omega)}\Vert v\Vert_{H^1_0(\Omega)}=M\Vert v\Vert_{H^1_0(\Omega)}.
  \end{aligned}$$

#### Nonhomogeneous Dirichlet condition

A problem with prescribed, nonzero boundary values can easily be
transformed into the case of zero boundary conditions. Suppose
$\partial \Omega$ is $C^1$ and $u\in H^1(\Omega)$, $f\in L^2(\Omega)$,
let $g$ be the trace of some $H^1$ function say $w$ in the trace sense,
see theorem ([8](#thm:trace){reference-type="ref"
reference="thm:trace"}). Then $u$ is a weak solution of
$$\label{eq:nonhomo}
  \begin{cases}
    -\Delta u&= f\textmd{ in }\Omega\\
    u&=g\textmd{ on }\partial \Omega,
  \end{cases}$$ which means $u=g$ on $\partial \Omega$ in the trace
sense and the identity
([\[eq:variation\]](#eq:variation){reference-type="ref"
reference="eq:variation"}) holds for all $v\in H^1_0(\Omega)$.

In order to make $u,v$ in same space, consider $\tilde u:= u-w$ belongs
to $H^1_0(\Omega)$, and consider $\tilde f\in H^{-1}(\Omega)$, the dual
space of $H_0^1$. Then the linear form in
([\[eq:variation\]](#eq:variation){reference-type="ref"
reference="eq:variation"}) be $$(\tilde f,v)=\int_\Omega fv - a(w,v),$$
and $\tilde u$ is a weak solution of the homogeneous boundary-value
problem $$\label{pro:non-dirichlet}
  \begin{cases}
    -\Delta\tilde u&= \tilde f\textmd{ in }\Omega\\
   \tilde u&=0\textmd{ on }\partial \Omega.
  \end{cases}$$

::: rem
**Remark 1**.
*$$H^1_0(\Omega)\subset L^2(\Omega)\subset H^{-1}(\Omega)$$*
:::

# Galerkin methods  {#sec:galerkin}

Consider the linear abstract variational problem: Find $u\in V$ such
that $$\forall v\in V, \quad a(u,v)=f(v),$$ wher the space $V$ and
bilinear form $a(\cdot, \cdot)$ and linear form $f$ are assumed to
satisfy the assumption of the Lax-Milgram Theorem.

Then the Galerkin method for approximating the solution of such a
problem sonsists in defining similar problems in finite-dimensional
subspaces of the space $V$.

::: defn
**Definition 2**. *Given finite-dimensional subspace $V_h\subset V$, we
associate the discrete problem:*

*Find $u_h\in V_h$ such that
$$\forall v_h\in V_h, \quad a(u_h,v_h)=f(v_h).$$ The unique solution
$u_h$ called a *discrete solution*.*
:::

**Three basic aspects of conforming finite element methods**

1.  (FEM 1) The first aspect is that a triangulation $\mathscr{T}_h$ is
    established over the set $\bar \Omega$, i.e., the set $\bar \Omega$
    is subdivided into a finite number of subsets $K$, called finite
    elements, in such a way that the following properties are satisfied:

    1.  $\bar \Omega = \cup_{K\in \mathscr{T}_h}K.$

    2.  For each $K\in \mathscr{T}_h$, the set $K$ is closed and the
        interior $\mathring{K}$ is nonempty.

    3.  For each distinct $K_1,K_2\in \mathscr{T}_h,$ one has
        $\mathring{K}_1\cap \mathring{K}_2=\emptyset.$

    4.  For each $K\in \mathscr{T}_h,$ the boundary $\partial K$ is
        Lipschitz-continuous.

2.  (FEM 2) The second basic aspect of the finite element method is that
    the spaces $P_K, K\in \mathscr{T}_h$, contain polynomials, or, at
    least contain functions which are \"close to\" polynomials.

    1.  it is the key to all convergence results as we shall see, and

    2.  it yields simple computations of the coefficients of the
        resulting linear system. (See below).

        Let $(w_k)_{k=1}^M$ be a basis in the space $V_h$. Then the
        solution $u_h=\sum_{k=1}^Mu_kw_k$ of problem $a(u_h,v_h)=f(v_h)$
        is such that the coefficients $u_k$ are solutions of the linear
        system
        $$\sum_{k=1}^Ma(w_k,w_l)u_k= f(w_l), \quad 1\leq l \leq M.$$
        whose matrix is always invertible, since the bilinear form,being
        assumed to be $V$-elliptic, is a fortiori $V_h$-elliptic.

3.  (FEM 3) There exists at least one \"canonical\" basis in the space
    $V_h$ whose corresponding basis functions have supports which are as
    \"small\" as possible.

## Finite elements {#sec:finite-ele}

There are lots of different ways to construct the triangulation and
spaces. Follow the main ideas above, we give a simple example of finite
elements here.

::: flushright
:::

Suppose, $\Omega\in \mathbb{R}^n$ is a polyhedral domain. Use
$n$-simplicial complex to give a triangulation ${\cal T}_h(\Omega)$.
Then a finite element in $\mathbb{R}^n$ is a triple $(K,P,\Sigma)$
where:

1.  $K$ is a $n$-simplex of $\mathbb{R}^n$,

2.  $P$ is a space of $d$-polynomials defined over the set $K$, with
    points $x_i,1\leq i\leq {n+d \choose d}$.

3.  $\Sigma$ is a finite set of nodal basis functions $\phi_i$,
    $1\leq i\leq N, N ={n+d \choose d}$, defined over the space $P$.

Therefore, for any function $v_h|_K\in P$ we have the representation
$$v_h|_K(x) = \sum _{i=1}^Nv_h(x_i)\phi_i(x).$$ Given a finite element
$(K,P,\Sigma)$ and given a function $v:K\rightarrow \mathbb{R}$ smooth,
define a project $\Pi$
$$\Pi: C^\infty(K)\rightarrow P(K)\quad  \Pi v= \sum_{i=1}^N\phi_i(v)v(x_i).$$

::: defn
**Definition 3**. *$$\begin{aligned}
      &h_K= \textmd{ diamenter of }K,\\
      &\rho_K=\textmd{ supremum of the diameters of the spheres inscribed in K,}\\
      &h= \max_{K\in \mathscr{T}_h} h_K.
    \end{aligned}$$*
:::

::: defn
**Definition 4**. *The associated finite element space $X_h$ of problem
([\[eq:nonhomo\]](#eq:nonhomo){reference-type="ref"
reference="eq:nonhomo"}) is the subspace of the product space
$\prod_{k\in \mathscr{T}_h}P_K$ defined by
$$S^d_h=\left\{v_h\in H^1(\Omega):\  v|_{K}\in P_K, v_h|_{\partial \Omega}=\Pi_hg
    \right\}$$ where $P_K$ is a space of $d$-th polynomials .*
:::

Though $v_h$ is a polynomial and thus $v_h\in H^1(K)$, it may has
problem on the common boundary $\Gamma$ of two simplex $P_1$ and $P_2$.
In order to make $v_h\in H^1(\Omega)$ we need $v_h$ be continuous on the
common boundaries. Since the boundary $\Gamma\in \mathbb{R}^{n-1}$, make
enough common points $x_i$ on $\Gamma$ we can make $v_h$ continuous on
$\Gamma.$ When $n=2$ and $P$ is a space of 2-polynomials, we only need
one more point, say mid-point to make $v_h$ continuous.

::: center
:::

# Analysis {#sec:analysis}

The variational equations $$\forall v\in V, \quad a(u,v) =f(v),$$
satisfy the assumptions of the Lax-Milgram Theorem.

With each finite element space $V_h$ is associated the discrete solution
$u_h$ which satisfies $$\forall v_h\in V_h,\quad a(u_h,v_h)=f(v_h).$$

## Convergence analysis {#sec:convergence}

::: defn
**Definition 5**. *If $\forall f$ of the variational problem, one has
$$\lim_{h\rightarrow 0}\Vert u-u_h\Vert =0,$$ Then we shall say that the
associated family of discrete problems is convergent.*
:::

::: defn
**Definition 6**. *If a family of finite elements satisfy there exists a
constant $\sigma$ such that $\forall h_K$, $h_K\leq \sigma \rho_K$ and
$h_K$ approach zero then we call it *a regular family of finite
elements.**
:::

## Stability {#sec:stability}

####  Nonhomogeneous dirichlet condition

Consider the homogeneous equation problem
([\[pro:non-dirichlet\]](#pro:non-dirichlet){reference-type="ref"
reference="pro:non-dirichlet"}) first.

Use Poincare's inequality for $\tilde u_h= (u-w)_h$ we have
$$\Vert \tilde u_h\Vert_{L^2(\Omega)}\leq C\Vert D\tilde u_h\Vert_{L^2(\Omega)}$$
then
$$\Vert \tilde u_h\Vert^2_{H^1(\Omega)}=\Vert \tilde u_h\Vert_{L^2(\Omega)}^2
  +\Vert D\tilde u_h\Vert_{L^2(\Omega)}^2
  \leq  C^2\Vert D\tilde u_h\Vert_{L^2(\Omega)}^2+ \Vert D\tilde u_h\Vert_{L^2(\Omega)}^2
  =(C^2+1) \Vert D\tilde u_h\Vert_{L^2(\Omega)}^2,$$ in addition,
$$\begin{aligned}
    \Vert D\tilde u_h\Vert_{L^2(\Omega)}^2
    &=a(\tilde u_h,\tilde u_h)
      =\tilde f(\tilde u_h)\\
    &=\int_{\Omega}f\cdot  \tilde u_h-a(w, \tilde u_h)\\
    &=\int_{\Omega}f\cdot  \tilde u_h-\int_{\Omega}Dw\cdot D\tilde u_h\\
    &\leq \Vert f\Vert_{L^2(\Omega)}\Vert \tilde u_h\Vert_{L^2(\Omega)}
      +\Vert Dw\Vert_{L^2(\Omega)}\Vert D\tilde u_h\Vert_{L^2(\Omega)}\\
    &\leq \Vert f\Vert_{L^2(\Omega)} \Vert \tilde u_h\Vert^2_{H^1(\Omega)}
      +\Vert Dw\Vert_{L^2(\Omega)} \Vert \tilde u_h\Vert^2_{H^1(\Omega)}
  \end{aligned}$$ hence, $$\Vert \tilde u_h\Vert_{H^1(\Omega)}\leq
   \Vert f\Vert_{L^2(\Omega)} +\Vert Dw\Vert_{L^2(\Omega)}$$ and the
discrete solution is bounded by the given condition so the method is
stable.

## Error equation

#### Estimate of the error

Use $\Pi_h u$ denote the interpolation of $u$ in $V_h^d$, with all
$v_h|_K$ is $d$-polynomial. Then the interpolation $\Pi_h u$ defined by
$u$ has interpolation error $$\Vert u-\Pi _hu\Vert _{H^1(\Omega)}
  =\left (\sum_{K\in \mathscr{T}_h}\Vert u-\Pi_h u\Vert^2_{H^1(\Omega)}\right)^{1/2}.$$
For all $v\in V_h\subset H_0^1(\Omega)$ we have
$$a(\tilde u,v)=\tilde f(v),\quad a(\tilde u_h,v)= \tilde f(v)$$ then
the error equation is $$a(\tilde u-\tilde u_h,v)=0.$$ Since
$\tilde u-\tilde u_h= \tilde u-\Pi_h \tilde u +\Pi_h \tilde u-\tilde u_h$,
we have $$\Vert \tilde u-\tilde u_h\Vert_{H^1(\Omega)}
  \leq \Vert \tilde u-\Pi_h \tilde u\Vert_{H^1(\Omega)}  +\Vert\Pi_h \tilde u-\tilde u_h\Vert_{H^1(\Omega)}$$
and $$a(\tilde u-\Pi_h \tilde u,v)+a(\Pi_h \tilde u-\tilde u_h,v)=0.$$
Let $v=\Pi_h \tilde u -\tilde u_h$, then
$$\beta \Vert v\Vert_{H^1(\Omega)}^2\leq a(v,v)$$ and
$$a(\tilde u-\Pi_h \tilde u,v)\leq
  \alpha \Vert v\Vert _{H^1(\Omega)}
  \Vert \tilde u-\Pi_h \tilde u\Vert _{H^1(\Omega)}$$ thus
$$\Vert \Pi_h \tilde u -\tilde u_h\Vert _{H^1(\Omega)}
  \leq C\Vert \tilde u-\Pi_h \tilde u\Vert _{H^1(\Omega)}$$ and
$$\Vert \tilde u-\tilde u_h\Vert_{H^1(\Omega)}\leq (C+1)\Vert \tilde u-\Pi_h \tilde u\Vert _{H^1(\Omega)}.$$
With all polynomials is $d$-order we have the interpolation error in
semi-norm that
$\vert \tilde u-\Pi_h \tilde u\vert _{0,\Omega}=O(h^{d+1})$ and
$\vert \tilde u-\Pi_h \tilde u\vert _{1,\Omega}=O(h^d)$, hence
$$\Vert \tilde u-\tilde u_h\Vert_{H^1(\Omega)}=O(h^d).$$

# Other examples {#sec:exam}

## Nonhomogeneous Neumann boundary condition {#sec:neumann}

Model problem: $$\label{eq:neumann}
  \begin{cases}
    -\Delta u+u&= f\quad \textmd{ in }\Omega\\
    \frac{\partial u}{\partial n}&=g \quad \textmd{ on }\partial \Gamma
  \end{cases}$$ with open bounded subset $\Omega\in \mathbb{R}^n$ and
$f\in L^2(\Omega), g\in L^2(\partial\Omega)$.

Multiply each side with trial function $v\in C^\infty(\Omega)$ and
integral in $\Omega$ we have
$$\int_{\Omega}-\Delta u v+uv=\int_\Omega fv$$ and then
$$\int_{\Omega} \nabla u\cdot \nabla v +uv= \int_{\partial \Omega}\frac{\partial u}{\partial n}v +\int_\Omega fv
  =\int_{\partial \Omega}gv +\int_\Omega fv,$$ thus we have bilinear
form $H^1(\Omega)\times H^1(\Omega)\rightarrow \mathbb{R}$ that
$$a(u,v)= \int_{\Omega} \nabla u\cdot \nabla v+uv,$$ and linear form
$H^1(\Omega)\rightarrow \mathbb{R}$ that
$$f(v)=\int_{\partial \Omega}gv +\int_\Omega fv.$$

#### Existence and uniqueness

Since
$$|a(u,v)|\leq \int_{\Omega} |\nabla u\cdot \nabla v|+\int_\Omega |uv|
  \leq |u|_{1,\Omega}|u|_{1,\Omega}+|u|_{0,\Omega}|u|_{0,\Omega}
  \leq \Vert u\Vert_{H^1_0(\Omega)}\Vert v\Vert_{H^1_0(\Omega)},$$ and
$$a(u,u)=\int_{\Omega} |\nabla u|^2+\int_\Omega u^2=\Vert u\Vert_{H^1(\Omega)}^2,$$
then $a(u,v)$ satisfies the hypotheses of the Lax-Milgram Theorem, so
that there exists a unique weak solution $u\in H^1(\Omega)$ of the
variation problem $$a(u,v)=f(v),\quad \forall v\in H^1(\Omega).$$

#### Stability {#stability}

Solve the weak form in finite-dimension to get discrete solution $u_h$
that $$a(u_h,v_h)=f(v_h)\quad \forall v_h\in V_h.$$

$$\begin{aligned}
    \Vert u_h\Vert_{H^1(\Omega)}^2
    &= \int_\Omega u_h\cdot u_h+\nabla u_h\cdot \nabla u_h\\
                                &=a(u_h,u_h)=f(u_h)\\
    &=\int_\Omega fu_h +\int_{\partial \Omega} gu_h\\
    &\leq \Vert g\Vert_{L^2(\partial \Omega)}\Vert u_h\Vert_{L^2(\Omega)}
      +\Vert f\Vert_{L^2(\Omega)}\Vert u_h\Vert_{L^2(\Omega)}\\
    &\leq ( \Vert g\Vert_{L^2(\partial \Omega)}+\Vert f\Vert_{L^2(\Omega)})
      \Vert u_h\Vert_{H^1(\Omega)}
  \end{aligned}$$ hence
$$\Vert u_h\Vert_{H^1(\Omega)}\leq  \Vert g\Vert_{L^2(\partial \Omega)}+\Vert f\Vert_{L^2(\Omega)}.$$

#### Convergence and error equation

Since $a(u,v)$ satisfy Lax-Milgram condition, the error equation can be
same as dirichlet boundary condition problems, say for all
$v\in V_h\subset H^1(\Omega)$ we have
$$a( u,v)= f(v),\quad a( u_h,v)=  f(v),$$ then $$a(u- u_h,v)=0.$$ Since
$u-  u_h=   u-\Pi_h   u +\Pi_h   u-  u_h$, we have
$$\Vert   u-  u_h\Vert_{H^1(\Omega)}
  \leq \Vert   u-\Pi_h   u\Vert_{H^1(\Omega)}  +\Vert\Pi_h   u-  u_h\Vert_{H^1(\Omega)}$$
and $$a(  u-\Pi_h   u,v)+a(\Pi_h   u-  u_h,v)=0.$$

Let $v=\Pi_h   u -  u_h$, then
$$\beta \Vert v\Vert_{H^1(\Omega)}\leq a(v,v)$$ and
$$a(  u-\Pi_h   u,v)\leq
  \alpha \Vert v\Vert _{H^1(\Omega)}
  \Vert   u-\Pi_h   u\Vert _{H^1(\Omega)}$$ thus
$$\Vert \Pi_h   u -  u_h\Vert _{H^1(\Omega)}
  \leq C\Vert   u-\Pi_h   u\Vert _{H^1(\Omega)}$$ and
$$\Vert   u-  u_h\Vert_{H^1(\Omega)}\leq (C+1)\Vert   u-\Pi_h   u\Vert _{H^1(\Omega)}.$$
With all polynomials is $d$-order we have the interpolation error in
semi-norm that $\vert   u-\Pi_h   u\vert _{0,\Omega}=O(h^{d+1})$ and
$\vert   u-\Pi_h   u\vert _{1,\Omega}=O(h^d)$, hence
$$\Vert   u-  u_h\Vert_{H^1(\Omega)}=O(h^d).$$

## Robin boundary condition {#sec:robin}

Model problem: $$\label{eq:neumann}
  \begin{cases}
    -\Delta u&= f\quad \textmd{ in }\Omega\\
    \frac{\partial u}{\partial n}+ \gamma u&=g \quad \textmd{ on }\partial \Gamma
  \end{cases}$$ with open bounded subset $\Omega\in \mathbb{R}^n$ and
$f\in L^2(\Omega), g\in L^2(\Omega)$, and $\gamma\neq 0$.

#### Weak solution {#weak-solution}

To define a weak solution we multiply the equation by $v\in H^1(\Omega)$
and integrate by parts to get $$\begin{aligned}
    \int_\Omega -\Delta uv&= \int_\Omega \nabla u\cdot \nabla v-\int_{\partial \Omega}\frac{\partial u}{\partial n}v\\
    &= \int_\Omega \nabla u\cdot \nabla v+\int_{\partial \Omega}\gamma uv
  \end{aligned}$$ then $u$ is a weak solution of
$$\int_\Omega \nabla u\cdot \nabla v+\gamma\int_{\partial \Omega} (Tu)(Tv)=\int_\Omega fv$$
where $u,v\in H^1(\Omega)$ and $T$ is trace operator, and that
$Tu=u|{\partial u}$. We define the linear form
$f:H^1(\Omega)\rightarrow \mathbb{R}$ similar as previous problem that
$$f(v)=\int_\Omega fv,$$ with Cauchy-Schwarz inequality, $f$ is
continuous. We define the bilinear form
$a(\cdot, \cdot):H^1(\Omega)\times H^1(\Omega)\rightarrow \mathbb{R}$ by
$$a(u,v)= \int_\Omega \nabla u\cdot \nabla v+\gamma\int_{\partial \Omega} (Tu)(Tv)\quad \forall u,v\in H^1(\Omega)$$

We prove that the bilinear form is continuous by the Cauchy-Schwarz
inequality and the Trace inequality form theorem
([8](#thm:trace){reference-type="ref" reference="thm:trace"}) that
$$\begin{aligned}
    |a(u,v)|&\leq \int_\Omega  \nabla u\cdot \nabla v+\gamma\int_{\partial \Omega} (Tu)(Tv)\\
            &\leq \Vert \nabla u\Vert_{L^2(\Omega)} \Vert \nabla v\Vert_{L^2(\Omega)}
              + \Vert Tu\Vert_{L^2(\partial\Omega)}\Vert Tv\Vert_{L^2(\partial\Omega)}\\
    &\leq \Vert \nabla u\Vert_{L^2(\Omega)} \Vert \nabla v\Vert_{L^2(\Omega)}
      + C^2\Vert u\Vert_{H^1(\Omega)}\Vert v\Vert_{H^1(\Omega)}\\
    & \leq(1+C^2) \Vert u\Vert_{H^1(\Omega)}\Vert v\Vert_{H^1(\Omega)}
  \end{aligned}$$ Then we need to consider the coercivity, i.e. we prove
the statement there exists constant $\beta$ such that
$$\beta\Vert u\Vert_{H^1(\Omega)}\leq a(u,u).$$

We argue by contradiction, that $a(\cdot, \cdot)$ is not coercive, then
the stated estimate false, there would exist for each integer
$k=1,\ldots,$ a function $u_k\in H^1(\Omega)$ satisfying $$\label{eq:1}
  \Vert u_k\Vert^2_{H^1(\Omega)}>k a(u_k,u_k).$$ We renormalize by
defining $$v_k=\frac{u_k}{\Vert u_k\Vert_{L^2(\Omega)}}.$$ Then
$$\Vert v_k\Vert_{L^p(\Omega)}=1,$$ and the equation
([\[eq:1\]](#eq:1){reference-type="ref" reference="eq:1"}) implies that
$$a(v_k,v_k)< \frac{1}{k}\quad (k=1,2,\ldots).$$ So $$\label{eq:2}
  \Vert \nabla u_k\Vert_{L^2(\Omega)}^2\leq a(u,u)< \frac{1}{k},\quad
  \Vert Tv_k\Vert_{L^2(\partial \Omega)}^2< \frac{1}{k}$$ and then
$\Vert v_k\Vert^2_{H^1(\Omega)}< 1+\frac{1}{k}.$ There exist a
subsequence $\{v_{k_j}\}_{j\geq 1}\subset \{v_k\}_{k\geq 1}$ and a
function $v\in L^2(\Omega)$ such that
$$v_{k_j}\rightarrow v\quad \textmd{ in }L^2(\Omega).$$ It follows that
$$\Vert v\Vert_{H^1(\Omega)}= 1.$$ On the other hand,
([\[eq:2\]](#eq:2){reference-type="ref" reference="eq:2"}) implies for
$i=1,2,\ldots,n$ and $\phi\in C^\infty_0(\Omega)$ that
$$\int_\Omega u\phi_{x_i} \textmd{d}x=
  \lim_{k_j\rightarrow\infty}\int_{\Omega} v_{k_j}\phi_{x_i}\textmd{d}x
  =- \lim_{k_j\rightarrow\infty}\int_{\Omega} v_{k_j,x_i}\phi\textmd{d}x=0$$
Consequently $v\in H^1(\Omega)$ with $\nabla v=0$. Thus $v$ is constant
on each component of $\Omega$. But
$$\Vert Tv\Vert_{L^2(\partial \Omega)}
  \leq \Vert T(v-v_k)\Vert_{L^2(\partial \Omega)}+\Vert Tv_k\Vert_{L^2(\partial \Omega)}
  \leq \Vert T\Vert \Vert v-v_k\Vert_{L^2(\partial \Omega)}+\sqrt{\frac{1}{k}}
  \rightarrow 0\quad \textmd{ as }k\rightarrow \infty.$$ which means
$v\equiv 0$ on $\partial \Omega$ and hence $v\equiv 0$ on $U$. This
contradicts to $\Vert v\Vert_{H^1}=1$. Therefore $a(\cdot,\cdot)$ is
coercive.

The Lax-Milgram theorem then gives a unique weak solution.

#### Stability {#stability-1}

Use the coercive of $a(\cdot, \cdot)$ we have
$$\beta\Vert u_h\Vert_{H^1(\Omega)}^2
  \leq a(u_h,u_h)= f(u_h)
  \leq \Vert f\Vert_{L^2(\Omega)}\Vert u_h\Vert_{L^2(\Omega)}
  \leq \Vert f\Vert_{L^2(\Omega)}\Vert u_h\Vert_{H^1(\Omega)}$$ and then
$$\Vert u_h\Vert_{H^1(\Omega)}\leq C \Vert f\Vert_{L^2(\Omega)}.$$

#### Error equation

Similar as previous question since $a(\cdot, \cdot)$ and $f$ satisfies
condition of Lax-Milgram Theorem.

# Sobolev spaces {#sec:sobolev}

::: {#thm:approx .thm}
**Theorem 7** (Global approximation by funcitons smooth up to the
boundary). *Assume $U$ is bounded and $\partial U$ is $C^1$. Suppose
$u\in W^{k,p}(U)$ for some $1\leq p <\infty$. Then there exist functions
$u_m\in C^\infty(\bar U)$ such that
$$u_m\rightarrow u\quad \textmd{ in }W^{k,p}(U).$$*
:::

::: {#thm:trace .thm}
**Theorem 8** (Trace Theorem). *Assume $U$ is bounded and $\partial U$
is bounded and $\partial U$ is $C^1$. Then there exists a bounded linear
operator $$T: W^{1,p}(U)\rightarrow L^p( \partial U)$$ such that*

1.  *$Tu=u|_{\partial U}$ if $u\in W^{1,p}(U)\cap C(\bar U)$*

2.  *$\Vert Tu\Vert _{L^p(\partial U)}\leq C\Vert u\Vert_{W^{1,p}(U)}$*

*for each $u\in W^{1,p}(U)$, with the constant $C$ depending only on $p$
and $U$.*
:::

::: {#eq:lax .thm}
**Theorem 9** (Lax-Milgram Theorem). *Let $V$ be a Hilbert space, let
$a(\cdot,\cdot):V\times V\rightarrow \mathbb{R}$ be a continuous
$V$-elliptic bilinear form, which means there exist constant
$\alpha,\beta>0$ such that
$$|a(u,v)|\leq \alpha \Vert u\Vert \Vert v\Vert,\quad (u,v\in V)$$ and
satisfy the coercivity condition
$$\beta\Vert u\Vert^2\leq a(u,u)\quad (u\in V).$$ and let
$f:V\rightarrow \mathbb{R}$ be a continuous linear form.*

*Then the abstract variational problem: Find an element $u$ such that
$$\label{eq:l}
    u\in V\quad \textmd{ and }\quad \forall v\in V,\quad  a(u,v)=f(v)$$
has one and only one solution.*
:::

::: thm
**Theorem 10**. *(Céa's lemma) There exists a constant $C$ independent
upon the subspace $V_h$ such that
$$\Vert u-u_h\Vert \leq C\inf_{v_h\in V_h}\Vert u-v_h\Vert.$$*
:::
