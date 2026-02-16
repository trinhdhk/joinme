// Some definitions presupposed by pandoc's typst output.
#let blockquote(body) = [
  #set text( size: 0.92em )
  #block(inset: (left: 1.5em, top: 0.2em, bottom: 0.2em))[#body]
]

#let horizontalrule = line(start: (25%,0%), end: (75%,0%))

#let endnote(num, contents) = [
  #stack(dir: ltr, spacing: 3pt, super[#num], contents)
]

#show terms: it => {
  it.children
    .map(child => [
      #strong[#child.term]
      #block(inset: (left: 1.5em, top: -0.4em))[#child.description]
      ])
    .join()
}

// Some quarto-specific definitions.

#show raw.where(block: true): set block(
    fill: luma(230),
    width: 100%,
    inset: 8pt,
    radius: 2pt
  )

#let block_with_new_content(old_block, new_content) = {
  let d = (:)
  let fields = old_block.fields()
  fields.remove("body")
  if fields.at("below", default: none) != none {
    // TODO: this is a hack because below is a "synthesized element"
    // according to the experts in the typst discord...
    fields.below = fields.below.abs
  }
  return block.with(..fields)(new_content)
}

#let empty(v) = {
  if type(v) == str {
    // two dollar signs here because we're technically inside
    // a Pandoc template :grimace:
    v.matches(regex("^\\s*$")).at(0, default: none) != none
  } else if type(v) == content {
    if v.at("text", default: none) != none {
      return empty(v.text)
    }
    for child in v.at("children", default: ()) {
      if not empty(child) {
        return false
      }
    }
    return true
  }

}

// Subfloats
// This is a technique that we adapted from https://github.com/tingerrr/subpar/
#let quartosubfloatcounter = counter("quartosubfloatcounter")

#let quarto_super(
  kind: str,
  caption: none,
  label: none,
  supplement: str,
  position: none,
  subrefnumbering: "1a",
  subcapnumbering: "(a)",
  body,
) = {
  context {
    let figcounter = counter(figure.where(kind: kind))
    let n-super = figcounter.get().first() + 1
    set figure.caption(position: position)
    [#figure(
      kind: kind,
      supplement: supplement,
      caption: caption,
      {
        show figure.where(kind: kind): set figure(numbering: _ => numbering(subrefnumbering, n-super, quartosubfloatcounter.get().first() + 1))
        show figure.where(kind: kind): set figure.caption(position: position)

        show figure: it => {
          let num = numbering(subcapnumbering, n-super, quartosubfloatcounter.get().first() + 1)
          show figure.caption: it => {
            num.slice(2) // I don't understand why the numbering contains output that it really shouldn't, but this fixes it shrug?
            [ ]
            it.body
          }

          quartosubfloatcounter.step()
          it
          counter(figure.where(kind: it.kind)).update(n => n - 1)
        }

        quartosubfloatcounter.update(0)
        body
      }
    )#label]
  }
}

// callout rendering
// this is a figure show rule because callouts are crossreferenceable
#show figure: it => {
  if type(it.kind) != str {
    return it
  }
  let kind_match = it.kind.matches(regex("^quarto-callout-(.*)")).at(0, default: none)
  if kind_match == none {
    return it
  }
  let kind = kind_match.captures.at(0, default: "other")
  kind = upper(kind.first()) + kind.slice(1)
  // now we pull apart the callout and reassemble it with the crossref name and counter

  // when we cleanup pandoc's emitted code to avoid spaces this will have to change
  let old_callout = it.body.children.at(1).body.children.at(1)
  let old_title_block = old_callout.body.children.at(0)
  let old_title = old_title_block.body.body.children.at(2)

  // TODO use custom separator if available
  let new_title = if empty(old_title) {
    [#kind #it.counter.display()]
  } else {
    [#kind #it.counter.display(): #old_title]
  }

  let new_title_block = block_with_new_content(
    old_title_block, 
    block_with_new_content(
      old_title_block.body, 
      old_title_block.body.body.children.at(0) +
      old_title_block.body.body.children.at(1) +
      new_title))

  block_with_new_content(old_callout,
    block(below: 0pt, new_title_block) +
    old_callout.body.children.at(1))
}

// 2023-10-09: #fa-icon("fa-info") is not working, so we'll eval "#fa-info()" instead
#let callout(body: [], title: "Callout", background_color: rgb("#dddddd"), icon: none, icon_color: black, body_background_color: white) = {
  block(
    breakable: false, 
    fill: background_color, 
    stroke: (paint: icon_color, thickness: 0.5pt, cap: "round"), 
    width: 100%, 
    radius: 2pt,
    block(
      inset: 1pt,
      width: 100%, 
      below: 0pt, 
      block(
        fill: background_color, 
        width: 100%, 
        inset: 8pt)[#text(icon_color, weight: 900)[#icon] #title]) +
      if(body != []){
        block(
          inset: 1pt, 
          width: 100%, 
          block(fill: body_background_color, width: 100%, inset: 8pt, body))
      }
    )
}



#let article(
  title: none,
  subtitle: none,
  authors: none,
  date: none,
  abstract: none,
  abstract-title: none,
  cols: 1,
  lang: "en",
  region: "US",
  font: "libertinus serif",
  fontsize: 11pt,
  title-size: 1.5em,
  subtitle-size: 1.25em,
  heading-family: "libertinus serif",
  heading-weight: "bold",
  heading-style: "normal",
  heading-color: black,
  heading-line-height: 0.65em,
  sectionnumbering: none,
  toc: false,
  toc_title: none,
  toc_depth: none,
  toc_indent: 1.5em,
  doc,
) = {
  set par(justify: true)
  set text(lang: lang,
           region: region,
           font: font,
           size: fontsize)
  set heading(numbering: sectionnumbering)
  if title != none {
    align(center)[#block(inset: 2em)[
      #set par(leading: heading-line-height)
      #if (heading-family != none or heading-weight != "bold" or heading-style != "normal"
           or heading-color != black) {
        set text(font: heading-family, weight: heading-weight, style: heading-style, fill: heading-color)
        text(size: title-size)[#title]
        if subtitle != none {
          parbreak()
          text(size: subtitle-size)[#subtitle]
        }
      } else {
        text(weight: "bold", size: title-size)[#title]
        if subtitle != none {
          parbreak()
          text(weight: "bold", size: subtitle-size)[#subtitle]
        }
      }
    ]]
  }

  if authors != none {
    let count = authors.len()
    let ncols = calc.min(count, 3)
    grid(
      columns: (1fr,) * ncols,
      row-gutter: 1.5em,
      ..authors.map(author =>
          align(center)[
            #author.name \
            #author.affiliation \
            #author.email
          ]
      )
    )
  }

  if date != none {
    align(center)[#block(inset: 1em)[
      #date
    ]]
  }

  if abstract != none {
    block(inset: 2em)[
    #text(weight: "semibold")[#abstract-title] #h(1em) #abstract
    ]
  }

  if toc {
    let title = if toc_title == none {
      auto
    } else {
      toc_title
    }
    block(above: 0em, below: 2em)[
    #outline(
      title: toc_title,
      depth: toc_depth,
      indent: toc_indent
    );
    ]
  }

  if cols == 1 {
    doc
  } else {
    columns(cols, doc)
  }
}

#set table(
  inset: 6pt,
  stroke: none
)

#set page(
  paper: "us-letter",
  margin: (x: 1.25in, y: 1.25in),
  numbering: "1",
)

#show: doc => article(
  title: [JoinME Stan Model Architecture],
  subtitle: [Comprehensive Technical Documentation],
  toc_title: [Table of contents],
  toc_depth: 3,
  cols: 1,
  doc,
)

#horizontalrule

= Overview
<overview>
The #strong[JoinME] package implements a sophisticated joint model for multivariate longitudinal and time-to-event data using Stan. This document provides an in-depth technical explanation of the Stan model architecture, data flow, and computational algorithms.

== Core Capabilities
<core-capabilities>
+ #strong[Multivariate Longitudinal Submodel];: Handles multiple markers measured repeatedly over time with flexible distribution families (Gaussian, Student-t, Bernoulli, Poisson, negative binomial, beta, ordinal, skew-normal, double exponential)
+ #strong[Survival Submodel];: Proportional hazards model with competing risks and flexible baseline hazard (B-spline)
+ #strong[Association Structure];: Links longitudinal and survival processes through:
  - Current value (CV) and current slope (CS) of longitudinal trajectories
  - Covariance structure of random effects
  - Flexible transformations of association features
+ #strong[Dynamic Prediction];: Conditional prediction of future longitudinal trajectories and survival probabilities given observed history

Additional implementation details that are central to the interpretation:

- Marker-weighted contributions are formed as weighted means across markers, using either fixed weights or shrinkage-estimated weights.
- The VCOV association depends on an id-level latent effect that drives a covariance regression, producing subject-specific Cholesky factors.
- Non-linear transforms are applied to CV, CS, and VCOV before they enter the hazard, using identity, functional opcodes, monotone I-splines, or piecewise linear rules.
- Competing risks are supported through event-type indexing; the overall survival in prediction integrates the sum of cause-specific hazards.

= High-Level Architecture
<high-level-architecture>
```{mermaid}
%%| fig-width: 100%
flowchart TD
    subgraph Data["Data Layer"]
        LongData[("Longitudinal Data<br/>N observations<br/>D markers<br/>n_id subjects")]
        EventData[("Survival Data<br/>n_id subjects<br/>K_event risks")]
        Covariates[("Covariates<br/>Fixed effects<br/>Hazard covariates<br/>Vcov covariates")]
    end
    
    subgraph Model["Model Layer"]
        direction TB
        FitModel["joinme_fit.stan<br/><b>Model Fitting</b>"]
        DynPred["joinme_dynpred.stan<br/><b>Dynamic Prediction</b>"]
    end
    
    subgraph Parameters["Parameter Space"]
        direction LR
        FixedEffects["β: Fixed Effects"]
        RandomEffects["u, v, w: Random Effects"]
        Hazard["γ, h0: Hazard Parameters"]
        Association["α: Association Weights"]
    end
    
    subgraph Likelihood["Likelihood Functions"]
        LongLik["Longitudinal Likelihood<br/>P(Y|RE, β)"]
        SurvLik[" Survival Likelihood<br/>P(T, δ|RE, γ, h0, α)"]
    end
    
   subgraph Helpers["Helper Functions"]
        PartialFit["partial_joinme()<br/>Subject-specific log-likelihood"]
        PartialDynPred["partial_draw()<br/>Draw-specific prediction"]
        Transform["Association Transformations<br/>apply_transform_*()"]
        CumHaz["Cumulative Hazard<br/>cumhaz()"]
        EtaFD["Finite Difference Slopes<br/>eta_fd()"]
    end
    
    LongData --> FitModel
    EventData --> FitModel
    Covariates --> FitModel
    
    FitModel --> Parameters
    Parameters --> DynPred
    LongData -.->|New Subject| DynPred
    
    FitModel --> PartialFit
    DynPred --> PartialDynPred
    
    PartialFit --> LongLik
    PartialFit --> SurvLik
    
    PartialDynPred --> LongLik
    PartialDynPred --> SurvLik
    
    LongLik --> Transform
    SurvLik --> Transform
    SurvLik --> CumHaz
    SurvLik --> EtaFD
    
    style FitModel fill:#e1f5ff
    style DynPred fill:#fff4e1
    style PartialFit fill:#ffe1f5
    style PartialDynPred fill:#f5e1ff
```

= Model Fitting: `joinme_fit.stan`
<model-fitting-joinme_fit.stan>
== Data Flow
<data-flow>
```{mermaid}
%%| fig-width: 100%
graph LR
    subgraph Input["Data Block"]
        LongObs["Longitudinal<br/>y_real, y_int<br/>X_obs, Z_*_obs"]
        SurvObs["Survival<br/>S_event, d_event<br/>W, Bs_*"]
        Design["Design Matrices<br/>X_gk_*, Z_*_gk_*<br/>Association info"]
    end
    
    subgraph TransData["Transformed Data"]
        IDRanges["id_start[i], id_end[i]<br/>Row ranges per subject"]
        Indices["Covariance indices<br/>idx_row_cov, idx_col_cov"]
    end
    
    subgraph Params["Parameters"]
        Beta["β: Fixed effects<br/>(scaled)"]
        RE["u_id: ID REs<br/>v_marker: Marker REs<br/>w_idscaled: Interaction REs"]
        Scales["τ: RE scales<br/>Lcorr: Correlations"]
        Gamma["γ: Hazard params<br/>h0: Baseline hazard"]
        Alpha["α: Association<br/>coefficients"]
    end
    
    subgraph Likelihood["Model Block"]
        ReduceSum["reduce_sum()<br/>Parallelizes over subjects"]
        PartialJoinME["partial_joinme(id_seq, start, end, ...)<br/>Likelihood slice"]
    end
    
    subgraph SubjectLik["Subject i Likelihood"]
        LongLoop["For each obs n ∈ [id_start[i], id_end[i]]:<br/>Compute η = Xβ + Zu + Zv + Zw<br/>Accumulate log p(y[n] | η, family)"]
        SurvCalc["Compute association features:<br/>CV, CS at GK nodes & event time<br/>Transform via functional ops<br/>Compute log h(t) & cumulative hazard"]
    end
    
    Input --> TransData
    TransData --> Likelihood
    Params --> Likelihood
    
    Likelihood --> ReduceSum
    ReduceSum -.->|Slice| PartialJoinME
    PartialJoinME --> SubjectLik
    SubjectLik --> LongLoop
    SubjectLik --> SurvCalc
    
    style ReduceSum fill:#ffe6e6
    style PartialJoinME fill:#fff0e6
```

== Algorithm Pseudocode
<algorithm-pseudocode>
=== Main Model Block
<main-model-block>
```plaintext
FOR each MCMC iteration:
  1. Sample parameters from priors
  2. Compute scaled fixed effects: β_scaled = β / tmax^{idx_time_*}
  3. Call reduce_sum(partial_joinme, id_seq, grainsize, ...)
     ├─ Splits subjects into chunks of size grainsize
     ├─ Processes each chunk in parallel (if threading enabled)
     └─ Returns total log-likelihood
  4. Add log-likelihood to target
END FOR
```

=== `partial_joinme()` Function
<partial_joinme-function>
This is the core likelihood evaluation per subject slice.

```plaintext
FUNCTION partial_joinme(id_seq, start, end, ... [all data & parameters]):
  lp = 0  // Accumulated log-likelihood for this slice
  
  FOR ii = start TO end:  // Process subjects in this slice
    i = id_seq[ii]  // Current subject index
    
    // -------------------
    // LONGITUDINAL LIKELIHOOD
    // -------------------
    FOR each observation n IN [id_start[i], id_end[i]]:
      d = marker[n]  // Marker type for this observation
      
      // 1. Compute linear predictor (location parameter)
      η_long = X_obs[n]'β_scaled             // Fixed effects
             + Z_id_obs[n]'u_id[i]           // ID-specific RE
             + Z_mk_obs[n]'v_marker[d]       // Marker-specific RE
             + Z_idm_obs[n]'w_idscaled[i,d]  // ID×Marker interaction RE
      
      // 2. Compute distributional parameters if needed
      //    (e.g., σ for Gaussian, ν for Student-t, φ for NegBinom, etc.)
      //    These may have their own regression structures (X_sigma, etc.)
      
      // 3. Accumulate log-density based on family
      SWITCH family_long[d]:
        CASE Gaussian:
          σ = exp(X_sigma[n]'β_sigma) OR sigma_marker[d] OR sigma_y
          lp += log N(y_real[n] | η_long, σ)
          
        CASE Student_t:
          σ = ... (as above)
          ν = 2 + exp(X_nu[n]'β_nu) OR nu_marker[d]
          lp += log t_ν(y_real[n] | η_long, σ)
          
        CASE Bernoulli:
          lp += log Bernoulli(y_int[n] | logit^{-1}(η_long))
          
        CASE Poisson:
          lp += log Poisson(y_int[n] | exp(η_long))
          
        CASE NegBinom:
          φ = exp(X_phi[n]'β_phi) OR phi_nb_marker[d]
          lp += log NegBinom(y_int[n] | exp(η_long), φ)
          
        ... [other families: beta, ordinal, skew_normal, etc.]
      END SWITCH
    END FOR
    
    // -------------------
    // SURVIVAL LIKELIHOOD
    // -------------------
    
    // 1. Compute association features at 15 Gauss-Kronrod (GK) quadrature nodes
    //    These nodes tile [0, S_event[i]] for numerical integration
    
    FOR j = 1 TO 15:  // For each GK node
      // a. Compute current value (CV) components
      cvm_now[j] = X_gk_now[i][j]'β_scaled + Z_id_gk_now[i][j]'u_id[i]  // Mean component
      cvk_now[j] = Z_mk_gk_now[i][j]'vbar + Z_idm_gk_now[i][j]'wbar_i[i]  // Marker component
      
      // b. Compute forward step for finite-difference slope
      cvm_fwd[j] = X_gk_fwd[i][j]'β_scaled + Z_id_gk_fwd[i][j]'u_id[i]
      cvk_fwd[j] = Z_mk_gk_fwd[i][j]'vbar + Z_idm_gk_fwd[i][j]'wbar_i[i]
    END FOR
    
    // 2. Compute current slopes (CS) via finite differences
    csm_raw = eta_fd(cvm_now, cvm_fwd, eps_fd) / tmax  // Time-rescaled slope
    csk_raw = eta_fd(cvk_now, cvk_fwd, eps_fd) / tmax

    // 2b. Marker-weighted means used for cvk and csk
    //     vbar and wbar_i are weighted means across markers using weights ω_d
    
    // 3. Compute variance-covariance association feature
    //    (logarithm of variances from Cholesky factor L_i)
    //    L_i is built by covariance regression with id-level latent effect u_L[i]
    //    lp_m = alpha_L[m] + beta_L[m]'Xcov[i] + lambda_L[m] * u_L[i]
    //    L_i[r,c] = log1p_exp(lp_m) if r == c, else lp_m
    vcov_raw = weighted_log_variances(L_i[i], a_vcov_var)
    
    // 4. Aggregate association features
    cv_total[j] = cvm_now[j] + cvk_now[j]  // Total current value
    cs_total[j] = csm_raw[j] + csk_raw[j]  // Total current slope
    
    // 5. Apply nonlinear transformations to association features
    //    Transformations can be: identity, log, exp, splines, functional opcodes, etc.
    cv_total_tf = apply_transform(cv_total, tf_mode_cv_tot, ops_cv, coeff_cv, ...)
    cv_mean_tf  = apply_transform(cvm_now, tf_mode_cv_mean, ...)
    cv_marker_tf = apply_transform(cvk_now, tf_mode_cv_marker, ...)
    cs_total_tf = apply_transform(cs_total, tf_mode_cs_tot, ...)
    cs_mean_tf  = apply_transform(csm_raw, tf_mode_cs_mean, ...)
    cs_marker_tf = apply_transform(csk_raw, tf_mode_cs_marker, ...)
    vcov_tf     = apply_transform(vcov_raw, tf_mode_vcov, ...)
    
    // 6. Combine transformed features into survival predictor at GK nodes
    FOR j = 1 TO 15:
      η_assoc_nodes[j] = α_cv_total * cv_total_tf[j]
                       + α_cv_mean * cv_mean_tf[j]
                       + α_cv_marker * cv_marker_tf[j]
                       + α_cs_total * cs_total_tf[j]
                       + α_cs_mean * cs_mean_tf[j]
                       + α_cs_marker * cs_marker_tf[j]
                       + vcov_tf[j]
    END FOR

    // 7. Competing risks handling
    //    For K_event causes, compute cause-specific hazards and sum in the
    //    cumulative hazard. The event-time contribution uses the observed cause.
    
    // 7. If subject had event (d_event[i] = 1), add log-hazard at event time
    IF d_event[i] == 1:
      k_ev = event_type[i]  // Which competing risk
      
      // Baseline hazard at event time
      log_h0_S = log_h0_intercept[k_ev] + Bs_event_c[i]'bs_gamma_c[k_ev]
      
      // Hazard covariate effects
      η_w_S = W[i]'γ_w[k_ev]
      
      // Association feature at event time (same steps as GK nodes)
      [Compute cv_*, cs_*, vcov_* at S_event[i]]
      [Transform and aggregate into η_assoc_S]
      
      // Add log-hazard contribution
      lp += log_h0_S + η_w_S + η_assoc_S
    END IF
    
    // 8. Compute cumulative hazard via Gauss-Kronrod quadrature
    //    Λ(t) = ∫₀ᵗ Σₖ hₖ(s) ds ≈ weighted sum over GK nodes
    
    FOR j = 1 TO 15:
      // Total hazard at node j across all K_event competing risks
      log_h_total[j] = log Σₖ exp(log_h0[k] + Bs_gk_c[i][j]'bs_gamma_c[k] + η_w + η_assoc_nodes[j])
    END FOR
    
    // Apply Gauss-Kronrod formula with weights w[j]
    Λ_i = cumhaz(S_event[i], log_h_total, zeros(15))
    
    // 9. Add survival contribution: log S(t) = -Λ(t)
    lp += -Λ_i
    
  END FOR  // End subject loop
  
  RETURN lp
END FUNCTION
```

=== Key Computational Steps Explained
<key-computational-steps-explained>
==== 1. #strong[Random Effects Structure]
<random-effects-structure>
Three-level random effects: $ u_(i d) [i] tilde.op upright("MVN") (0 \, Sigma_(i d)) $ $ v_(m a r k e r) [d] tilde.op upright("MVN") (0 \, Sigma_(m a r k e r)) $ $ w_(i d times m) [i \, d] tilde.op upright("MVN") (0 \, Sigma_(i d m)) $

Covariance matrices parameterized via Cholesky: $ Sigma_(i d) = upright("diag") (tau_(i d)) dot.op L_(i d)^(c o r r) (L_(i d)^(c o r r))' dot.op upright("diag") (tau_(i d)) $ $ Sigma_(m a r k e r) = upright("diag") (tau_(m k)) dot.op L_(m k)^(c o r r) (L_(m k)^(c o r r))' dot.op upright("diag") (tau_(m k)) $ $ Sigma_(i d m) = L_i L_i' $ where $L_i$ is a subject-specific Cholesky factor driven by covariance regression.

$macron(w) [i]$ is the weighted average of $w$ across markers (for "total" associations), and $macron(v)$ is the weighted average of $v$ across markers.

==== 2. #strong[Association Mechanisms]
<association-mechanisms>
#strong[Current Value (CV)];: - $C V_(t o t a l) = upright("Fixed") + u_(i d) + v_(m k) + w_(i d times m)$ - $C V_(m e a n) = upright("Fixed") + u_(i d)$ (subject-average trajectory) - $C V_(m a r k e r) = v_(m k) + w_(i d times m)$ (marker-specific deviations)

#strong[Current Slope (CS)];: - Computed via finite differences: $ C S (t) approx frac(C V (t + epsilon) - C V (t), epsilon) $

#strong[Covariance (VCOV)];: - Weighted log-variances from Cholesky factor $L_i$. - Captures uncertainty in subject-specific predictions.

==== 3. #strong[Gauss-Kronrod Quadrature]
<gauss-kronrod-quadrature>
To integrate hazard over $[0 \, T]$: $ Lambda (T) = integral_0^T lambda (s) d s $

Use 15-point Gauss-Kronrod rule: $ Lambda (T) approx T / 2 sum_(j = 1)^15 omega_j lambda (T dot.op x_j) $ where $x_j in [0 \, 1]$ are GK nodes and $omega_j$ are weights. Design matrices are pre-computed at these nodes.

==== 4. #strong[Flexible Transformations]
<flexible-transformations>
Association features can be transformed before entering survival predictor. Transformation modes: - #strong[0];: Identity ($f (x) = x$) - #strong[1];: Functional opcodes (log, exp, sqrt, compositions) - #strong[2];: I-spline basis (for monotonic smooth transforms) - #strong[3];: Piecewise linear

Example: If CV exhibits nonlinear association with hazard: $ eta_(a s s o c) + = alpha dot.op f (C V) $ where $f ()$ is learned via splines.

= Dynamic Prediction: `joinme_dynpred.stan`
<dynamic-prediction-joinme_dynpred.stan>
== Conceptual Overview
<conceptual-overview>
Given: - Posterior samples {θ^(k)} from fitted model - New subject with observed history y\_obs, times t\_obs - Conditioning time T\_cond (usually max(t\_obs))

Goal: - Predict future longitudinal values y\_future | y\_obs, θ - Predict survival probability S(t | T\_cond, y\_obs, θ)

Competing risks note: When K\_event \> 1, the prediction returns the overall survival probability that integrates the sum of cause-specific hazards. Cause-specific survival curves are not returned directly by the default dynpred outputs.

Strategy: 1. For each posterior draw k: a. Sample subject-specific random effects u, v, w | y\_obs, θ^(k) b. Compute longitudinal and survival predictions given u, v, w 2. Average predictions across draws to get marginal estimates

== Data Flow
<data-flow-1>
```{mermaid}
%%| fig-width: 100%
graph TB
    subgraph Input["Inputs"]
        PostDraws["Posterior Draws<br/>β^(k), τ^(k), γ^(k), α^(k)<br/>k = 1...n_draws"]
        NewSub["New Subject Data<br/>y_obs, X_obs, Z_obs<br/>Observed history"]
        PredGrid["Prediction Grids<br/>t_future (longitudinal)<br/>t_surv (survival)"]
    end
    
    subgraph Model["Model Block"]
        ReduceDraw["reduce_sum()<br/>Parallelizes over draws"]
        PartialDraw["partial_draw(draw_ids, start, end, ...)<br/>Draw-specific likelihood"]
    end
    
    subgraph DrawK["For Draw k"]
        SampleRE["Sample Random Effects<br/>u ~ N(0,I), v ~ N(0,I), w_lat ~ N(0,I)<br/>Transform: u ← L_u·z_u, etc."]
        BuildCov["Build Subject Cholesky L_i<br/>from covariance regression"]
        LongCond["Longitudinal Conditioning<br/>log p(y_obs | u,v,w, θ^(k))"]
        SurvCond["Survival Conditioning<br/>log S(T_cond | u,v,w, θ^(k))"]
    end
    
    subgraph GenQuant["Generated Quantities"]
        PredLong["Predict y_future<br/>at each t_future"]
        PredSurv["Predict S(t | T_cond)<br/>at each t_surv"]
        Summarize["Aggregate over draws<br/>Mean, Median, CrI"]
    end
    
    PostDraws --> Model
    NewSub --> Model
    Model --> ReduceDraw
    ReduceDraw -.->|Slice| PartialDraw
    PartialDraw --> DrawK
    DrawK --> SampleRE
    SampleRE --> BuildCov
    BuildCov --> LongCond
    BuildCov --> SurvCond
    LongCond -.-> GenQuant
    SurvCond -.-> GenQuant
    GenQuant --> PredLong
    GenQuant --> PredSurv
    GenQuant --> Summarize
    
    style PartialDraw fill:#fff0e6
    style SampleRE fill:#e6f7ff
```

== Algorithm Pseudocode
<algorithm-pseudocode-1>
=== `partial_draw()` Function
<partial_draw-function>
This function operates on a slice of posterior draws.

```plaintext
FUNCTION partial_draw(draw_ids_slice, start, end, ... [data & posterior draws]):
  lp = 0  // Accumulated log-density for this draw slice
  
  FOR i = start TO end:  // Process draws in slice
    k = draw_ids_slice[i]  // Current draw index
    
    // -------------------
    // 1. SAMPLE RANDOM EFFECTS FROM STANDARD NORMAL
    // -------------------
    z_u[k]      ~ N(0, I_{R_id})            // ID-level standard normal
    z_v[k,d]    ~ N(0, I_{R_mk})  ∀d        // Marker-level standard normal
    z_w_lat[k,d] ~ N(0, I_{Q_idm}) ∀d       // Latent ID×marker standard normal
    z_L[k]      ~ N(0, 1)                   // Latent covariance regression
    
    lp += log p(z_u, z_v, z_w_lat, z_L)  // Standard normal densities
    
    // -------------------
    // 2. TRANSFORM TO ACTUAL RANDOM EFFECTS
    // -------------------
    // Build Cholesky factors from draw k parameters
    L_u = diag(τ_id^(k)) * Lcorr_id^(k)
    L_v = diag(τ_marker^(k)) * Lcorr_marker^(k)
    L_w = diag(τ_marker_id^(k)) * Lcorr_marker_id^(k)
    
    // Transform standard normals to correlated REs
    u_id = L_u * z_u[k]
    
    FOR each marker d:
      v_marker[d] = L_v * z_v[k, d]
    END FOR
    
    // For subject×marker interaction, allow cross-marker correlation
    FOR each marker d:
      cross = B_cross^(k) * v_marker[d]  (if cross-correlation enabled)
      z_w[d] = cross + L_w * z_w_lat[k, d]
    END FOR
    
    // -------------------
    // 3. BUILD SUBJECT-SPECIFIC COVARIANCE VIA REGRESSION
    // -------------------
    // Latent random intercept for covariance
    u_Lk = τ_vcov_reg^(k) * z_L[k]
    
    // Build Cholesky factor L_i for this subject
    L_i = matrix(Q_idm × Q_idm, 0)
    FOR each unique entry (r,c) in lower triangle:
      // Regression on covariance covariates
      lin_pred = α_vcov_reg^(k)[r,c]
               + vec_cov_vcov' * β_vcov_reg^(k)[r,c]
               + λ_vcov_reg^(k)[r,c] * u_Lk
      
      IF r == c:  // Diagonal (log scale with softplus)
        L_i[r,c] = log(1 + exp(lin_pred))
      ELSE:       // Off-diagonal (unrestricted)
        L_i[r,c] = lin_pred
      END IF
    END FOR
    
    // Transform latent REs to actual subject×marker REs
    FOR each marker d:
      w_idscaled[d] = L_i * z_w[d]
    END FOR

    // 3b. Marker-weighted means for association summaries
    //     weights ω_d are draw-specific when marker weights are estimated
    vbar = Σ_d ω_d v_marker[d]
    wbar = Σ_d ω_d w_idscaled[d]
    
    // -------------------
    // 4. LONGITUDINAL CONDITIONING: log p(y_obs | u, v, w, θ^(k))
    // -------------------
    FOR each observed longitudinal obs n:
      d = idx_marker_obs[n]
      
      // Linear predictor
      η = mat_fixed_obs[n]' * β_fixed^(k)
        + mat_id_obs[n]' * u_id
        + mat_marker_obs[n]' * v_marker[d]
        + mat_marker_id_obs[n]' * w_idscaled[d]
      
      // Likelihood based on family (same logic as joinme_fit)
      SWITCH family_long[d]:
        CASE Gaussian:
          σ = ... (from β_sigma^(k) or σ_marker^(k))
          lp += log N(y_real[n] | η, σ)
        ... [other families]
      END SWITCH
    END FOR
    
    // -------------------
    // 5. COMPUTE ASSOCIATION FEATURES AT GK NODES (conditional on T_cond)
    // -------------------
    // Same logic as joinme_fit, but at conditioning time T_cond
    
    // Weighted averages of REs across markers
    marker_weights_k = marker_weights_draws[k, :]
    vbar = weighted_average(v_marker, marker_weights_k)
    wbar_i = L_i * weighted_average(z_w, marker_weights_k)
    
    FOR j = 1 TO 15:  // GK nodes in [0, T_cond]
      cvm[j] = mat_fixed_gk_cond[j]' * β_fixed^(k) + mat_id_gk_cond[j]' * u_id
      cvm_f[j] = mat_fixed_gk_cond_fwd[j]' * β_fixed^(k) + mat_id_gk_cond_fwd[j]' * u_id
      
      cvk[j] = mat_marker_gk_cond[j]' * vbar + mat_marker_id_gk_cond[j]' * wbar_i
      cvk_f[j] = mat_marker_gk_cond_fwd[j]' * vbar + mat_marker_id_gk_cond_fwd[j]' * wbar_i
    END FOR
    
    // Compute slopes, aggregate, transform
    csm_raw = eta_fd(cvm, cvm_f, eps_finite_diff)
    csk_raw = eta_fd(cvk, cvk_f, eps_finite_diff)
    vcov_raw = weighted_log_variances(L_i, coeff_assoc_vcov_var^(k))
    
    cv_total = cvm + cvk
    cs_total = csm_raw + csk_raw
    
    // Apply transformations
    cv_total_tf = apply_transform(cv_total, ...)
    ... [other transforms]
    
    // Aggregate association
    η_assoc_nodes = α_cv_total^(k) * cv_total_tf
                  + α_cv_mean^(k) * cv_mean_tf
                  + ... [other terms]
    
    // -------------------
    // 6. SURVIVAL CONDITIONING: log S(T_cond | u, v, w, θ^(k))
    // -------------------
    FOR j = 1 TO 15:
      // Total hazard at node j
      log_h_total[j] = log Σₖ exp(log_h0^(k)[k] + mat_basis_gk_cond[j]' * bs_gamma_c^(k)[k]
                                  + vec_cov_hazard' * γ_hazard^(k)[k]
                                  + η_assoc_nodes[j])
    END FOR
    
    // Cumulative hazard up to T_cond
    Λ_cond = cumhaz(time_condition, log_h_total, zeros(15))
    
    // Add survival contribution
    lp += -Λ_cond
    
  END FOR  // End draw loop
  
  RETURN lp
END FUNCTION
```

=== Generated Quantities
<generated-quantities>
Once random effects are sampled conditional on observed data, the `generated quantities` block:

+ #strong[Longitudinal Predictions];:

  ```plaintext
  FOR each future time t_fut:
    FOR each marker d:
      η_fut[d, t_fut] = X_fut[t_fut]'β + Z_id_fut[t_fut]'u + Z_mk_fut[t_fut]'v[d] + Z_idm_fut[t_fut]'w[d]

      // Expectation (epred scale)
      y_pred_mean[d, t_fut] = inv_link(η_fut, family[d])

      // Predictive draw (with residual noise)
      y_pred[d, t_fut] ~ distribution(η_fut, family[d])
    END FOR
  END FOR
  ```

+ #strong[Survival Predictions];:

  ```plaintext
  FOR each future time t_surv > T_cond:
    // Compute cumulative hazard from T_cond to t_surv
    Λ_future = cumhaz_from_to(T_cond, t_surv, log_h_total, ...)

    // Conditional survival
    S(t_surv | T_cond) = exp(-Λ_future)
  END FOR
  ```

+ #strong[Aggregation];:

  ```plaintext
  Across all draws k = 1...n_draws:
    mean_y_pred[d, t_fut] = mean_k(y_pred_mean[k, d, t_fut])
    quantiles_y_pred[d, t_fut] = quantile_k(y_pred_mean[k, d, t_fut], [0.025, 0.5, 0.975])

    mean_S(t | T_cond) = mean_k(S[k, t | T_cond])
    quantiles_S(t | T_cond) = quantile_k(S[k, t | T_cond], [0.025, 0.5, 0.975])
  ```

= Helper Functions Documentation
<helper-functions-documentation>
== `eta_fd(now, fwd, eps)`: Finite Difference Slope
<eta_fdnow-fwd-eps-finite-difference-slope>
#strong[Purpose];: Approximate the derivative (slope) of a trajectory using forward finite differences.

#strong[Algorithm Setup];: - `now`: vector of values at time $t$ - `fwd`: vector of values at time $t + epsilon$ - `eps`: step size $epsilon$

#strong[Computation];: $ upright("slope") [j] = frac(upright("fwd") [j] - upright("now") [j], epsilon) $

#strong[Usage in Joint Model];: - Computes "current slope" (CS) association features. - Applied at Gauss-Kronrod nodes and event times. - Time rescaling by $1 \/ T_(m a x)$ converts to original time scale.

== `cumhaz(S, log_h0_nodes, eta_nodes)`: Cumulative Hazard via Gauss-Kronrod
<cumhazs-log_h0_nodes-eta_nodes-cumulative-hazard-via-gauss-kronrod>
#strong[Purpose];: Numerically integrate hazard function over $[0 \, S]$ using 15-point Gauss-Kronrod quadrature.

#strong[Algorithm Setup];: - $S$: upper limit of integration (scaled event/censoring time) - `log_h0_nodes`: log baseline hazard at 15 GK nodes - `eta_nodes`: hazard predictor at 15 GK nodes

#strong[Output];: $ Lambda (S) = integral_0^S h (t) d t $

#strong[Computation];: Using weights $omega = [0.0229 \, 0.0631 \, dots.h]$: 1. Compute log-hazard contributions: $z_j = log h_0 (x_j) + eta (x_j)$ 2. Numerically stable log-sum-exp accumulation: $ upright("acc") = sum_(j = 1)^15 omega_j dot.op exp (z_j - max (z)) $ 3. Scale result: $ Lambda (S) = S / 2 dot.op exp (max (z)) dot.op upright("acc") $

== `eta_vcov_varonly_weighted_const(K, L_i, a_vcov_var)`: Variance Association
<eta_vcov_varonly_weighted_constk-l_i-a_vcov_var-variance-association>
#strong[Purpose];: Compute weighted log-variances from a Cholesky decomposition, used as an association feature.

#strong[Algorithm Setup];: - $K$: number of evaluations - $L_i$: $Q times Q$ Cholesky factor ($Sigma_i = L_i L_i'$) - $a_(v c o v)$: $Q$-vector of association weights

#strong[Computation];: $ s = sum_(q = 1)^Q a_(v c o v) [q] dot.op log (sum_(k = 1)^q L_i [q \, k]^2) $ Returns a constant vector of length $K$.

== `apply_transform_*()`: Flexible Nonlinear Transformations
<apply_transform_-flexible-nonlinear-transformations>
#strong[Purpose];: Apply user-specified transformations to association features before entering survival predictor.

#strong[Supported Modes];:

+ #strong[Identity] (`mode = 0`): $f (x) = x$

+ #strong[Functional Opcodes] (`mode = 1`):

  - Stack-based evaluator for compositions: $f (x) = log (1 + exp (x))$
  - Operations: log, exp, sqrt, inv\_logit, logit, square, reciprocal, add/mult const, power, min/max, abs, trig functions.

+ #strong[I-Spline Basis] (`mode = 2`):

  - Monotonic smooth transformation via integrated B-splines.
  - Enforces $f' (x) gt.eq 0$.

+ #strong[Piecewise Linear] (`mode = 3`):

  - Linear interpolation between knot points.

#strong[Algorithm] (functional opcodes): Evaluates reverse polish notation Opcodes on a stack. For example `PUSH_X, PUSH_1, ADD, LOG` computes $log (x + 1)$.

#strong[Vectorized Version];: - For vectors, apply element-wise (same opcodes, broadcast constants) - Efficient for evaluating at multiple GK nodes simultaneously

= Computational Considerations
<computational-considerations>
== Threading and Parallelization
<threading-and-parallelization>
=== `reduce_sum()` in `joinme_fit.stan`
<reduce_sum-in-joinme_fit.stan>
```plaintext
target += reduce_sum(partial_joinme, id_seq, grainsize, ...)

grainsize: Number of subjects per parallel chunk
- grainsize = n_id → serial execution (single chunk)
- grainsize = ceil(n_id / n_threads) → parallel execution

Threading model:
- Each chunk processed independently
- No shared mutable state → thread-safe
- Results summed after all chunks complete
```

#strong[Trade-offs];: - #strong[Small grainsize];: More parallelism, higher overhead - #strong[Large grainsize];: Less parallelism, lower overhead - #strong[Optimal];: Depends on n\_id, model complexity, hardware

=== Threading in `joinme_dynpred.stan`
<threading-in-joinme_dynpred.stan>
```plaintext
target += reduce_sum(partial_draw, draw_ids, grainsize, ...)

Parallelizes over posterior draws instead of subjects
- Each draw's likelihood computed independently
- Particularly beneficial when n_draws is large (e.g., 200-1000)
```

== Memory and Efficiency
<memory-and-efficiency>
=== Design Matrix Pre-computation
<design-matrix-pre-computation>
All association design matrices (X\_gk\_#emph[, Z\_];\_gk\_\*, etc.) are pre-computed in R and passed as data. This avoids: - Recomputing spline bases at every MCMC iteration - Building design matrices inside Stan loops - Memory allocations during sampling

#strong[Trade-off];: Larger data object, but much faster sampling.

=== Sparse Representations
<sparse-representations>
Random effects designs (Z\_\*) are typically sparse. Stan's matrix operations leverage this implicitly, but explicit sparse matrices are not used due to compatibility constraints with `reduce_sum`.

=== Caching Cholesky Factors
<caching-cholesky-factors>
Correlation matrices parameterized via Cholesky factors (Lcorr\_\*) to: - Ensure positive-definiteness - Avoid costly eigendecompositions - Enable efficient sampling via standard normal draws

== Numerical Stability
<numerical-stability>
=== Log-Sum-Exp for Competing Risks
<log-sum-exp-for-competing-risks>
$ log (sum_(k = 1)^K exp (h_k)) = max (h) + log (sum_(k = 1)^K exp (h_k - max (h))) $

Prevents overflow/underflow when hazards have vastly different magnitudes.

=== Softplus for Positive Parameters
<softplus-for-positive-parameters>
$ sigma = log (1 + exp (eta)) $

Alternative to $exp (eta)$ when avoiding extreme values (always $> 0$, smooth, stable).

=== Time Rescaling
<time-rescaling>
All times scaled to $[0 \, 1]$ via division by $t_(m a x)$: - Numerical stability (avoids large numbers in exponentials) - Coefficients interpretable on standardized scale - Fixed effects rescaled post-hoc for original scale

= Common Pitfalls and Debugging
<common-pitfalls-and-debugging>
== Divergent Transitions
<divergent-transitions>
#strong[Symptom];: Stan warnings about divergent transitions after warmup.

#strong[Causes] in JoinME context: 1. #strong[Stiff survival likelihood];: Hazard changes rapidly, creating challenging geometry 2. #strong[Weakly identified random effects];: Some REs not well-constrained by data 3. #strong[Prior-data conflict];: Priors incompatible with likelihood

#strong[Remedies];: - Increase `adapt_delta` (e.g., 0.95 or 0.99) → smaller step sizes - Reparameterize: Use non-centered parameterization for REs if not already - Tighten priors on poorly identified parameters - Simplify model: Remove weak associations or unnecessary REs

== Non-Finite Log-Likelihood
<non-finite-log-likelihood>
#strong[Symptom];: `log_prob = -inf` or `nan`.

#strong[Causes];: 1. #strong[Distributional parameter out of bounds];: e.g., σ ≤ 0, ν ≤ 0 for Student-t 2. #strong[Incorrect data encoding];: Integer outcomes for continuous family, or vice versa 3. #strong[Transformation domain violation];: log of negative, sqrt of negative

#strong[Debugging];: - Check `family_long` indices match data types - Verify all distributional regressions produce valid ranges (use link functions) - Add `print()` statements in Stan to trace parameter values

== Slow Sampling
<slow-sampling>
#strong[Symptom];: \< 10 iterations/second, or hours for 1000 iterations.

#strong[Causes];: 1. #strong[Large n\_id or N];: More subjects/observations = more computation 2. #strong[Complex associations];: Many transformations, high-dimensional REs 3. #strong[Inefficient grainsize];: Too small → overhead dominates

#strong[Remedies];: - Reduce `n_draws` for testing - Simplify model: Disable unused associations - Tune `grainsize`: Try ceil(n\_id / (4 \* threads)) - Use threading: `threads_per_chain = 2` or `4`

== Prediction Failures
<prediction-failures>
#strong[Symptom];: `no chains finished successfully` in `predict()`.

#strong[Causes];: 1. #strong[Threaded model not compiled with threading];: Exe missing `stan_threads` 2. #strong[Insufficient warmup for fixed\_param sampler];: `iter_warmup` too small 3. #strong[Incompatible data dimensions];: Prediction grid misaligned

#strong[Remedies];: - Clear Stan cache, recompile dynpred models - Increase `iter_warmup` to 100+ - Validate `times` vector is within valid range

= Advanced Topics
<advanced-topics>
== Covariance Regression
<covariance-regression>
Subject-specific Cholesky factor $L_i$ is modeled via regression. For each unique entry $(r \, c)$ of $L_i$:

$ eta_(r \, c) = alpha [r \, c] + X_(c o v) [i] beta [r \, c] + lambda [r \, c] u_L [i] $

If $r = c$ (Diagonal): $ L_i [r \, c] = upright("softplus") (eta_(r \, c)) $ Else (Off-diagonal): $ L_i [r \, c] = eta_(r \, c) $

#strong[Enables];: - Covariance structure varying with subject covariates. - Individual-specific uncertainty in longitudinal predictions ($u_L$). - Richer association with survival (via vcov feature).

#strong[Shrinkage];: - Normal priors on $beta [r \, c]$ with scale 0.3 - Random intercept $u_L$ $arrow.r$ subject-level variation.

== Marker Weighting
<marker-weighting>
Association features use #strong[weighted averages] of marker-specific REs:

$ macron(v) = frac(1, D) sum_d omega_d v_(m a r k e r) [d] $ $ macron(w)_i = L_i dot.op frac(1, D) sum_d omega_d z_w [d] $

#strong[Weights] $omega_d$: - Fixed: User-specified base weights (e.g., equal weights). -- Estimated: Signed additive perturbation around base weights with RMS stabilization, scaled by $tau_(m a r k e r \_ w e i g h t s)$. -- Reflects relative importance of each marker for survival after shrinkage while allowing negative and positive contributions (Park & Casella 2008).

#strong[Interpretation];: - High $omega_d$ $arrow.r$ marker $d$ strongly associated with event.

== Competing Risks
<competing-risks>
Hazard for subject $i$ is sum over $K_(e v e n t)$ risk types:

$ lambda (t \| i) = sum_(k = 1)^(K_(e v e n t)) lambda_k (t \| i) $ $ lambda_k (t \| i) = h_(0 k) (t) exp (W [i] ' gamma_k + eta_(a s s o c) (t \, i)) $

#strong[Event Type];: - If $d_(e v e n t) [i] = 1$, observed $e v e n t \_ t y p e [i] in { 1 \, dots.h \, K_(e v e n t) }$. - Log-likelihood includes $h_(e v e n t \_ t y p e [i]) (S_(e v e n t) [i])$. - Cumulative hazard includes all $K_(e v e n t)$ hazards.

#strong[Use Cases];: - Death vs.~dropout - Multiple disease endpoints - Recurrent events (with time-varying baseline)

== Distributional Regression
<distributional-regression>
Parameters like $sigma \, nu \, phi.alt$ can have their own regression structures:

$ log (sigma [n]) = X_sigma [n] beta_sigma + Z_sigma [n] b_sigma $

#strong[Allows];: - Heteroscedasticity (residual variance varies with covariates). - Marker-specific variance. - Random effects in variance (group-level variation).

#strong[Families Supporting Distributional Regression];: - Gaussian: σ (scale) - Student-t: σ, ν (scale, df) - Negative binomial: φ (dispersion) - Beta: φ (precision) - Skew-normal: σ, α (scale, skewness) - Skew-double-exponential: σ, τ (scale, skewness)

= Practical Workflow
<practical-workflow>
== 1. Model Fitting
<model-fitting>
```r
fit <- joinme(
  formulaLong = y ~ time + x1 + (1 + time || id) + (1 || marker),
  dataLong = long_data,
  formulaEvent = Surv(time, event) ~ x2,
  dataEvent = event_data,
  assoc = c("cv_total", "cs_total"),  # Association features
  family = "gaussian",  # Or auto-detect
  basehaz = "bs",  # B-spline baseline hazard
  n_knots = 5,
  control = list(
    chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.95,
    threads_per_chain = 2  # Enable threading
  )
)
```

== 2. Inspect Convergence
<inspect-convergence>
```r
summary(fit)       # Parameter estimates with R-hat, ESS
plot(fit)          # Traceplots
pp_check(fit)      # Posterior predictive checks
```

== 3. Dynamic Prediction
<dynamic-prediction>
```r
pred <- predict(
  fit,
  newdataLong = new_long,
  newdataEvent = new_event,
  Tstart = 5,  # Condition on data up to t=5
  times = seq(5, 10, length.out = 50),  # Predict t ∈ [5, 10]
  n_samples = 200
)

plot(pred, type = "survival")
plot(pred, type = "longitudinal", marker = "CD4")
```

= Priors and Regularization Strategy
<priors-and-regularization-strategy>
This section provides the definitive technical specification and rationale for the prior ecosystem within `joinme`. The prior choices strike a balance between #strong[weakly informative defaults] (allowing data to drive inference) and #strong[principled regularization] (preventing pathological behavior in complex joint likelihood surfaces).

== 1. Longitudinal Submodel
<longitudinal-submodel>
=== Fixed Effects ($beta$)
<fixed-effects-beta>
- #strong[Prior];: Student-$t (nu = 6 \, mu = 0 \, sigma = 2.5)$ (scaled).
- #strong[Rationale];: The Student-$t$ distribution with moderate degrees of freedom ($nu = 6$) provides fatter tails than a Normal distribution. This robustness allows for occasional large coefficients without overly penalizing them, while the scale prevents the sampler from exploring implausibly large regions of parameter space during warm-up.
- #strong[Support];: User-configurable scales allow adaptation to the magnitude of covariates.

=== Random Effect Scales ($tau_(i d) \, tau_(m k)$)
<random-effect-scales-tau_id-tau_mk>
- #strong[Prior];: Half-Normal($0 \, sigma$).
- #strong[Rationale];: The half-normal prior provides regularization toward zero (parsimony) while strictly enforcing positivity. Unlike the Half-Cauchy (often recommended in older literature), the Half-Normal has lighter tails, which reduces the risk of the "funnel" pathology where the sampler gets stuck in regions of extremely high variance where the data provide little information. This is crucial for convergence in high-dimensional random effects structures.

=== Random Effect Correlations ($Omega$)
<random-effect-correlations-omega>
- #strong[Prior];: LKJ-Correlation($eta$).
- #strong[Rationale];: The LKJ distribution over correlation matrices allows explicit control over the prior volume near the identity matrix.
  - $eta = 1$: Uniform over valid correlation matrices.
  - $eta > 1$ (Default $eta = 2$): Peaked at the identity matrix. This mildly regularizes correlations toward zero (independence), stabilizing estimation when clusters (subjects) are few or observation counts are low.
  - This avoids the degeneracy often seen with Inverse-Wishart priors on covariance matrices.

=== Distributional Parameters
<distributional-parameters>
Family-specific auxiliary parameters require constraints matching their support: - #strong[Scales ($sigma_y$)];: Exponential($lambda$). The maximum entropy prior for positive continuous variables with a known mean. - #strong[Degrees of Freedom ($nu$)];: Gamma($alpha \, beta$). Prior constraints $nu > 2$, soft-bounded to avoid infinite variance regimes while allowing normality ($nu arrow.r oo$). - #strong[Skew/Quantile Parameters];: Beta($alpha \, beta$). Natural conjugate for bounded $(0 \, 1)$ parameters, defaulting to Uniform(0,1) or centered symmetric Beta(2,2).

== 2. Survival Submodel
<survival-submodel>
=== Baseline Hazard ($h_0 (t)$)
<baseline-hazard-h_0t>
- #strong[Coefficients ($gamma_(b s)$)];: Normal($0 \, sigma_(b s)$).
- #strong[Smoothness Penalty];: When B-splines are used, `joinme` imposes a #strong[random-walk prior] (second-order difference penalty) on the spline coefficients. $ gamma_k - 2 gamma_(k - 1) + gamma_(k - 2) tilde.op cal(N) (0 \, tau_(s m o o t h)) $
- #strong[Rationale];: This effectively smooths the baseline hazard, preventing overfitting to small clusters of events (spikes) while allowing flexibility in the hazard shape over time. This is critical for robust extrapolation in dynamic prediction.

=== Hazard Covariates ($gamma_w$)
<hazard-covariates-gamma_w>
- #strong[Prior];: Student-$t (nu = 6 \, mu = 0 \, sigma = 2.5)$.
- #strong[Rationale];: Consistent with fixed effects; robust regularization.

== 3. Association and Link Functions
<association-and-link-functions>
The association module is the most complex component and requires hierarchical shrinkage to handle high-dimensional marker systems ($D gt.double 1$).

=== Association Coefficients ($alpha$)
<association-coefficients-alpha>
The total association is often a sum of a mean term, a marker-specific deviation, and interactions. - #strong[Base Effect ($alpha_(m e a n)$)];: Normal($0 \, 1$). Weakly informative. - #strong[Marker-Specific Deviations ($alpha_(m a r k e r)$)];: - #strong[Normal Shrinkage];: $alpha_d tilde.op cal(N) (0 \, tau_(s h r i n k))$. Equivalent to Ridge regression ($L_2$ penalty). Encourages all markers to share strength. - #strong[Laplace Shrinkage];: $alpha_d tilde.op upright("Laplace") (0 \, tau_(s h r i n k))$. Equivalent to Lasso regression ($L_1$ penalty). Encourages sparsity, effectively selecting relevant markers that drive the risk. - #strong[Rationale];: In joint models with many biomarkers, it is a priori unlikely that all markers contribute equally to the risk. Shrinkage priors automatically attenuate noise markers, improving out-of-sample prediction performance.

=== Marker Weights ($omega_d$)
<marker-weights-omega_d>
When estimating a composite score (e.g., weighted sum of current values), weights are modeled as signed perturbations around base weights and then stabilized to unit RMS. -- #strong[Prior];: Normal or Laplace shrinkage on perturbations (controlled by `shrinkage`), scaled by $tau_(m a r k e r \_ w e i g h t s)$ (Park & Casella 2008). - #strong[Rationale];: Allows protective and harmful marker effects while controlling non-identifiable global weight scale.

=== Nonlinear Transformations
<nonlinear-transformations>
When `type = "ispline_penalized"`, the association function $f (x)$ is modeled via Monotone I-splines. - #strong[Prior];: Second-order difference penalty strictly on the I-spline coefficients. - #strong[Constraint];: Coefficients must be positive to ensure monotonicity. - #strong[Rationale];: Nonlinear associations (e.g., risk increases only after a threshold) are powerful but prone to overfitting in sparse data regions. The penalty forces the function toward linearity in the absence of strong evidence for curvature.

=== Covariance Regression (VCOV)
<covariance-regression-vcov>
- #strong[Latent Effect ($u_L$)];: Normal(0, $tau_L$).
- #strong[Coefficients];: Normal priors on regression terms linking baseline covariates to the log-volatility.
- #strong[Rationale];: This structure allows the variability of the biomarkers (and their correlation) to itself be a predictor of survival (e.g., highly labile blood pressure predicts stroke), regularized by the hierarchical structure of the latent variable.

= References
<references>
#block[
] <refs>

#horizontalrule

#strong[Document Version];: 1.0 \
#strong[Last Updated];: #strong[?meta:date] \
#strong[Package];: JoinME \
#strong[Stan Version];: ≥ 2.32
