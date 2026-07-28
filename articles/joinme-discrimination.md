# Discrimination analysis: Concordance and Time-dependent AUC

## 1 Purpose and scope

This vignette gives a transparent account of the discrimination measures
implemented for a fitted `JoiNMeFit` object. It is intended to provide
enough detail for a statistical methods section, independent
reproduction of the estimators, and critical assessment of their
assumptions.

The two public functions answer different scientific questions:

- `concordance()` asks whether the fitted **conditional survival curves
  order observed failure times correctly over the available follow-up**.
  It returns one follow-up-wide index.
- [`auc()`](https://trinhdhk.github.io/joinme/reference/auc.md) asks
  whether the fitted **conditional event risk at a specified horizon
  separates cases from controls at a chosen landmark**. It returns one
  cumulative/dynamic area under the receiver operating characteristic
  curve for every requested landmark–horizon combination.

These are measures of **discrimination**, not calibration. A model may
order subjects correctly while systematically overestimating or
underestimating their absolute risks. Conversely, well-calibrated group
averages do not imply that individual subjects are ordered well.

This distinction is the reason that `concordance()` does not require
`time_horizon` or `Dt`, whereas
[`auc()`](https://trinhdhk.github.io/joinme/reference/auc.md) does.
Supplying a horizon changes the scientific question from follow-up-wide
concordance to horizon-specific classification.

### 1.1 One prediction pathway, two estimands

Both measures begin with the same fitted joint model and the same
prospective restriction of each longitudinal history. They separate only
after posterior mean conditional survival has been obtained.

Code

``` default
flowchart LR
    A["Fitted joint model<br/>and evaluation data"] -->
    B["One terminal outcome<br/>per subject"]
    B --> C["Choose landmark and retain<br/>history observed by that time"]
    C --> D["Predict posterior mean<br/>conditional survival"]
    D --> E["Concordance:<br/>compare curves at each<br/>observable event time"]
    D --> F["AUC:<br/>compare risks at one<br/>specified horizon"]
```

The shared prediction pathway is important: a difference between the two
reported numbers arises from their comparison sets and time summaries,
not from fitting two different survival models.

## 2 Statistical setting and notation

### 2.1 Observed event data

For subject i=1,\ldots,n, let

- xT_i^\ast be the latent event time;
- C_i be the censoring time;
- T_i=\min(T_i^\ast,C_i) be the observed event or censoring time;
- \delta_i=1 indicate an observed event and \delta_i=0 indicate
  censoring;
- K_i be the observed event cause when competing events are present;
- \mathcal H_i(s) be the longitudinal history observed through time s;
  and
- \mathbf X_i be the event-process covariates.

The fitted event formula may use either a right-censored
`Surv(time, status)` response or counting-process
`Surv(start, stop, status)` data. Counting-process data can contain
several rows for one subject. Discrimination nevertheless concerns one
terminal outcome per subject. JoiNMe therefore orders the intervals and
retains the terminal row before forming any case, control, or
comparable-pair indicator. Left- and interval-censored outcomes are not
silently reinterpreted; they are rejected because they require different
comparison rules.

### 2.2 Conditional prediction origin

Let T\_{0i} denote the time at which prediction is made for subject i.
Only longitudinal observations at or before T\_{0i} are allowed to
inform the prediction. Define residual follow-up by

R_i = T_i-T\_{0i}.

The conditional survival function on this residual time scale is

S_i(u\mid T\_{0i}) = \Pr\\\left( T_i^\ast\>T\_{0i}+u \mid
T_i^\ast\>T\_{0i}, \mathcal H_i(T\_{0i}), \mathbf X_i, \mathcal
D\_{\mathrm{fit}} \right), \qquad u\geq 0,

where \mathcal D\_{\mathrm{fit}} denotes the data used to obtain the
posterior distribution. By construction, S_i(0\mid T\_{0i})=1.

JoiNMe calculates survival over posterior draws and uses the posterior
mean ordinate

\widehat S_i(u\mid T\_{0i}) = \operatorname{E}\\\left\[ S_i(u\mid
T\_{0i},\boldsymbol\theta) \mid \mathcal D\_{\mathrm{fit}} \right\]

for discrimination. The present functions therefore report a plug-in
discrimination estimate from posterior mean curves. They do not report
the posterior distribution, standard error, or confidence interval of
the discrimination measure.

### 2.3 How the joint model enters a survival prediction

For cause k, the event model may be written schematically as

\lambda\_{ik}(t) = \lambda\_{0k}(t) \exp\\\left\\ \mathbf X_i(t)^\mathsf
T\boldsymbol\gamma_k + \mathcal A\_{ik}\\\left( t,\mathcal H_i(t)
\right) \right\\,

where \lambda\_{0k}(t) is the fitted baseline hazard and \mathcal
A\_{ik} is the complete longitudinal–event association contribution. The
resulting event-free survival has the form

S_i(t\mid s) = \exp\\\left\[ -\sum_k\int_s^t \lambda\_{ik}(v)\\dv
\right\], \qquad t\>s,

with the conditioning adjustment that survival equals one at s.

Consequently, the discrimination functions do not reconstruct a
simplified risk score. They call
[`predict.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/predict.JoiNMeFit.md)
and inherit:

- the fitted baseline-hazard basis and its original centring constants;
- event covariates and time-varying event-process rows;
- posterior random effects informed by the longitudinal history;
- marker weights;
- current-value, current-slope, functional, correlation, and covariance
  association terms; and
- monotone-spline and ordered piecewise-linear association transforms.

For an ordered piecewise-linear association, the ordered ordinates
generated by the simplex and cumulative sum enter the hazard at every
quadrature and prediction time. No separate approximation or
unconstrained coefficient is introduced by `concordance()` or
[`auc()`](https://trinhdhk.github.io/joinme/reference/auc.md).

## 3 Follow-up-wide survival-curve concordance

### 3.1 The estimand

Suppose subject i experiences the cause of interest at residual time
R_i, and subject j is known to remain event-free beyond that residual
time. Their observable order is

R_i\<R_j.

Following the survival-curve principle proposed by [Antolini, Boracchi,
and Biganzoli (2005)](https://doi.org/10.1002/sim.2427), JoiNMe
evaluates **both** conditional survival curves at the earlier observed
event time R_i. The pair is:

\begin{aligned} \text{concordant} \quad& \text{if}\quad \widehat
S_i(R_i\mid T\_{0i}) \< \widehat S_j(R_i\mid T\_{0j}),\\
\text{discordant} \quad& \text{if}\quad \widehat S_i(R_i\mid T\_{0i}) \>
\widehat S_j(R_i\mid T\_{0j}),\\ \text{a prediction tie} \quad&
\text{if}\quad \widehat S_i(R_i\mid T\_{0i}) = \widehat S_j(R_i\mid
T\_{0j}). \end{aligned}

Lower survival means greater predicted risk. Thus a concordant pair
assigns greater risk at R_i to the subject who actually experiences the
event at R_i.

Let \mathcal P be the set of observable ordered pairs and let a_i be the
event-time weight attached to comparisons generated by event i. Define

c\_{ij} = \mathbb I\\\left\\ \widehat S_i(R_i)\<\widehat S_j(R_i)
\right\\ + \frac{1}{2} \mathbb I\\\left\\ \widehat S_i(R_i)=\widehat
S_j(R_i) \right\\.

The reported concordance is

\widehat C = \frac{ \displaystyle \sum\_{(i,j)\in\mathcal P}a_i c\_{ij}
}{ \displaystyle \sum\_{(i,j)\in\mathcal P}a_i }.

A value of 0.5 represents chance ordering under continuous predictions,
a value approaching 1 represents increasingly correct ordering, and a
value below 0.5 indicates systematic reversal. These descriptions
concern ranking only and do not establish calibration or clinical
usefulness.

### 3.2 Comparable pairs and censoring

Censoring does not imply long survival. It only limits which temporal
orderings are observed. JoiNMe follows the risk-set convention described
for
[`survival::concordance()`](https://stat.ethz.ch/R-manual/R-devel/library/survival/html/concordance.html).
For an event of interest at residual time r:

| Outcome for comparator j | Comparable? | Reason |
|----|---:|----|
| Event or censoring after r | Yes | Subject j is known to survive beyond the event time |
| Censoring exactly at r | Yes | Under the event-before-censor convention, j is known to be event-free through r |
| Censoring before r | No | The order between the event and j’s latent event time is unknown |
| Event of interest exactly at r | No | The two failure times are tied and have no observed strict order |
| Competing event before r | No | It is treated as cause-specific censoring before the comparison time |
| Competing event exactly at r | Yes | It follows the same event-before-censor convention as other censoring at r |

This pairwise rule uses a prematurely censored subject whenever its
observed history is sufficient for a particular earlier comparison, and
excludes it after its censoring time. It neither declares the subject a
late survivor nor imputes an unobserved failure order.

The interpretation relies on the censoring assumptions appropriate to
the chosen weighting method. Informative loss to follow-up that remains
unexplained by the fitted data can bias discrimination estimates. No
comparison rule alone can identify the ordering of two latent event
times after observation has ceased.

### 3.3 Event-time weighting

The `type_weights` argument selects an event-time weighting schedule.
JoiNMe asks `survival::concordance(..., ranks = TRUE)` for its
documented event-time weights and divides each event weight by the
number at risk at that event time. This produces the pair weight a_i
used in the survival-curve comparisons.

| `type_weights` | Interpretation in JoiNMe |
|----|----|
| `"none"` or `"n"` | Harrell-style pair weighting; every observable pair has weight one |
| `"S"` | Weight events using the estimated survival schedule |
| `"S/G"` | Weight events using the survival and censoring schedules |
| `"n/G2"` | Uno-style inverse-censoring event weighting supplied by `survival` |
| `"I"` | Equal event-time weighting supplied by `survival` |

For `"n"`, the event-time weight returned by `survival` equals the
risk-set size, so division by that size gives unit weight to every pair.
The other options change the relative contribution of event times.

The `"n/G2"` option reuses the event-time weighting associated with
Uno’s concordance in `survival`; the subject scores in JoiNMe remain the
time-dependent survival-curve comparisons defined above. It should
therefore be described as an Uno-style censoring-weighted Antolini
comparison, not as a claim that the implementation is identical in every
respect to a conventional static-score Uno estimator.

### 3.4 Choice of prediction origin

By default, `concordance(fit)` chooses

T\_{0i} = \max\\t\<T_i: Y_i(t)\text{ was observed}\\.

The strict inequality prevents a measurement recorded at the outcome
time from informing a prediction assessed against that same outcome.
Subjects without a valid pre-outcome longitudinal history cannot
contribute.

This default answers a residual-prognosis question: **given each
subject’s most recent pre-outcome information, does the model order
their remaining follow-up?** Because T\_{0i} may differ between
subjects, residual times rather than absolute study times are compared.

A common landmark gives a different and often more directly comparable
estimand:

Code

``` r

concordance(fit, time_start = 2)
```

This asks whether predictions made at study time 2 order subsequent
event times among subjects still under observation after time 2. The
landmark may also be a named vector or the name of a subject-constant
column in the event data:

Code

``` r

subject_landmark <- c("1" = 1.5, "2" = 2.0, "3" = 1.0)

concordance(
  fit,
  newdataLong = validation_long,
  newdataEvent = validation_event,
  time_start = subject_landmark
)

concordance(
  fit,
  newdataLong = validation_long,
  newdataEvent = validation_event,
  time_start = "clinical_landmark"
)
```

Subject-specific origins are appropriate only when the scientific target
is remaining time from those origins. They should not be presented as
concordance on a common absolute study-time scale.

### 3.5 Cause-specific analysis and competing events

The `cause` argument identifies the failure type that generates
comparable event pairs. Other event causes are treated as censoring at
their observed times. They may serve as comparators before their
competing event, but not after it.

The current [`predict()`](https://rdrr.io/r/stats/predict.html) result
used as the score is the fitted event-free `Survival`, which reflects
the complete fitted event process. Consequently, with competing risks,
the pair construction is cause-specific but its ranking quantity is
overall event-free survival. This is transparent and reproducible, but
it is not the same estimand as a cause-specific cumulative-incidence
concordance. A scientific question about cumulative incidence in the
presence of competing risks requires a dedicated competing-risk
discrimination estimator.

## 4 Concordance: implementation from call to result

The following sequence describes every public and internal function used
by the concordance calculation. Internal functions begin with a full
stop and are not part of the stable user interface.

### 4.1 Step 1: validate the request and resolve the data

[`concordance.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/concordance.JoiNMeFit.md):

1.  verifies that `object` is a fitted `JoiNMeFit`;
2.  checks that `cause` is a single positive integer;
3.  asks `.concordance_time_weight()` to validate and translate
    `type_weights`;
4.  calls `.resolve_train_data()`; and
5.  passes the aligned data to `.concordance_survival_curves()`.

`.resolve_train_data()` uses the stored fitting data only when both
`newdataLong` and `newdataEvent` are omitted. If either new data set is
provided, both are required. This prevents an accidental mixture of
validation longitudinal histories with training event outcomes, or the
converse.

`.concordance_time_weight()` maps `"none"` to the `survival` spelling
`"n"` and rejects any unknown weighting rule before prediction begins.

### 4.2 Step 2: construct one terminal outcome per subject

`.subject_event_outcomes()`:

1.  finds the `Surv()` response in the fitted event formula;
2.  verifies that the subject identifier is present and non-missing;
3.  accepts right-censored and counting-process responses;
4.  orders counting-process intervals by subject, stop time, and start
    time;
5.  retains the final interval for each subject; and
6.  decodes the event indicator and event cause into `event_status` and
    `event_type`.

The ordering makes the reduction invariant to the row order supplied by
the user. Retaining one row avoids counting a subject once for every
time-dependent covariate interval.

### 4.3 Step 3: determine and validate landmarks

When `time_start` is omitted, `.last_preoutcome_measurement()` finds the
final longitudinal measurement strictly before each observed event or
censoring time. A numerical tolerance prevents a value that differs only
through floating-point representation from being treated as earlier.

When `time_start` is supplied, `.subject_time_map()` converts it to one
finite, named value per subject. It accepts:

- one scalar shared by all subjects;
- a named vector matched by subject identifier;
- an unnamed vector already in subject order; or
- a character value naming a subject-constant event-data column.

For counting-process event data, `.constant_finite_value()` verifies
that a named landmark column is finite and constant within subject.

Subjects are eligible only when their observed outcome occurs strictly
after their landmark. The cause-specific event indicator is then

d_i=\mathbb I(\delta_i=1,\\ K_i=\texttt{cause}).

### 4.4 Step 4: restrict the longitudinal history

`.concordance_survival_curves()` retains only measurements satisfying

t\_{i\ell}\leq T\_{0i}.

This is an essential prospective-prediction condition. Measurements
after the landmark could improve estimation of a subject’s random
effects, but they were not available when the prediction was notionally
made and would introduce future information. Subjects without any usable
history after this restriction are removed.

### 4.5 Step 5: build prediction times

For each subject, `.concordance_prediction_grid()` constructs an
absolute-time grid containing:

- residual time zero, expressed as the subject’s landmark;
- 50 regularly spaced points between the landmark and observed follow-up
  end;
- every distinct cause-specific residual event time at which that
  subject can be a comparator; and
- the subject’s observed follow-up end.

Including every required event time makes the concordance ordinate exact
with respect to the prediction grid. The regular points retain a stable
curve for other prediction summaries and protect against sparse
interpolation.

### 4.6 Step 6: obtain conditional survival predictions

`.concordance_survival_curves()` calls:

``` r

predict(
  object,
  process = "event",
  times = subject_specific_time_grids,
  time_start = subject_specific_landmarks,
  control = list(n_samples = n_samples),
  seed = seed
)
```

[`predict.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/predict.JoiNMeFit.md)
evaluates the fitted hazard and its integral for the selected posterior
draws. `n_samples` controls the number of posterior draws used to
estimate each posterior mean survival curve; `seed` makes draw
subsampling reproducible. Additional prediction arguments supplied
through `...` are passed on unchanged.

### 4.7 Step 7: align curves on residual event times

`.concordance_survival_matrix()` constructs a matrix whose rows are
distinct cause-specific residual event times and whose columns are
subjects. For every subject it:

1.  subtracts that subject’s landmark from the absolute prediction
    times;
2.  orders the residual grid;
3.  linearly interpolates at the required event times; and
4.  returns `NA` outside the subject’s predicted support.

The necessary event times were inserted explicitly in Step 5, so
interpolation normally returns a stored ordinate. It mainly guards
against harmless floating-point changes during time rescaling.
Prohibiting extrapolation ensures that a subject is not compared beyond
their observed support.

### 4.8 Step 8: obtain event weights

`.concordance_event_weights()` constructs a right-censored residual-time
response and calls `survival::concordance(..., ranks = TRUE)` with an
arbitrary fixed working score. The score is used only because the
`survival` interface requires one; it cannot affect the risk-set sizes
or event-time weights.

The function aligns the returned rank table to event rows, including the
single-event case in which R may drop a matrix dimension. It divides
each event-time weight by the number of subjects whose residual
follow-up reaches that time. Non-event rows receive weight zero.

### 4.9 Step 9: visit every observable pair once

`.concordance_from_survival_curves()` visits each cause-specific event
and selects comparators who:

- have residual follow-up greater than the event time; or
- are censored exactly at the event time; and
- have a finite survival ordinate at that event time.

It subtracts the event subject’s survival from every comparator’s
survival, counts positive, negative, and zero differences, and
multiplies the counts by the event’s pair weight. Each strict temporal
ordering is visited once. The final estimate gives half credit to a
predicted tie.

The implementation stores an event-time-by-subject matrix rather than a
posterior-draw-by-time-by-subject array. If there are E distinct event
times, its principal storage is O(En). Pair comparison is O(n^2) in the
worst case, which is unavoidable for an exact all-pairs empirical
concordance without further structure.

### 4.10 Step 10: return an auditable summary

The result contains:

| Column | Meaning |
|----|----|
| `concordance` | ( \text{concordant}+0.5 \times \text{tied} )/ n\_{\text{pairs}} |
| `concordant` | Weighted number of correctly ordered pairs |
| `discordant` | Weighted number of incorrectly ordered pairs |
| `tied` | Weighted number of survival-probability ties |
| `n_pairs` | Total comparison weight |
| `n_subjects` | Number of subjects retained after landmark and history checks |
| `n_events` | Number of observed events of the requested cause |

With non-unit event weighting, the first four pair-count columns are
weighted amounts and need not be integers.

## 5 Cumulative/dynamic time-dependent AUC

### 5.1 The estimand

Choose a common landmark s and horizon t\>s. For a subject known to be
event-free at s, define posterior mean conditional risk

\widehat r_i(s,t) = 1-\widehat S_i(t-s\mid s).

JoiNMe uses the cumulative/dynamic case and control definitions:

D_i(s,t) = \begin{cases} 1, & T_i\in(s,t\],\\ \delta_i=1,\\
K_i=\texttt{cause}, \\ 0, & T_i\>t, \\ \text{unknown}, &
\text{otherwise}. \end{cases}

Thus:

- a **case** experiences the cause of interest after the landmark and by
  the horizon;
- a **control** is observed to remain event-free beyond the horizon;
- a subject censored in (s,t\] has unknown case/control status and is
  omitted;
- a competing event in (s,t\] is treated as cause-specific censoring and
  is omitted; and
- a subject whose outcome occurred at or before s is outside the dynamic
  risk set.

For n_1 cases and n_0 controls, the empirical cumulative/dynamic AUC is

\widehat{\operatorname{AUC}}(s,t) = \frac{1}{n_1n_0}
\sum\_{i:D_i=1}\sum\_{j:D_j=0} \left\[ \mathbb I\\\widehat r_i\>\widehat
r_j\\ + \frac{1}{2}\mathbb I\\\widehat r_i=\widehat r_j\\ \right\].

It is the probability that a randomly selected observed case receives
greater predicted horizon risk than a randomly selected observed
control, with half credit for a risk tie.

### 5.2 Censoring interpretation

The current AUC is an unweighted empirical complete-case estimator. It
does **not** apply inverse-probability-of-censoring weights.
Early-censored subjects are omitted because their case/control state at
the horizon is unobserved.

This makes the calculation and its effective sample transparent, but it
can change the target population when censoring is substantial. It is
most defensible when censoring before the horizon is limited and
plausibly independent of the latent event status, conditionally on the
analysis design. For substantial or informative censoring, an IPCW or
model-based time-dependent AUC should be pre-specified instead. The
returned case, control, and pair counts should always be reported.

### 5.3 Efficient calculation by rank sum

Construct the combined vector of case risks followed by control risks
and let \mathcal R_i be its average rank, so tied values receive their
mean rank. The number of control risks below the cases, with half credit
for ties, is

\sum\_{i:D_i=1}\mathcal R_i-\frac{n_1(n_1+1)}{2}.

Therefore,

\widehat{\operatorname{AUC}}(s,t) = \frac{ \displaystyle
\sum\_{i:D_i=1}\mathcal R_i-n_1(n_1+1)/2 }{ n_1n_0 }.

This Mann–Whitney identity is numerically equivalent to enumerating all
n_1n_0 case–control comparisons, but requires only a combined vector and
its ranks rather than a full pairwise matrix.

## 6 AUC: implementation from call to result

### 6.1 Step 1: dispatch and argument validation

[`auc()`](https://trinhdhk.github.io/joinme/reference/auc.md) is the S3
generic. `AUC` is an alias retained for users who prefer the capitalised
spelling.
[`auc.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/auc.md):

1.  validates the fitted object, landmark values, event cause, and
    horizons;
2.  resolves training or new data through `.resolve_train_data()`;
3.  determines a horizon for every landmark; and
4.  processes each landmark–horizon pair in turn.

When `time_horizon` is omitted:

- `Dt` gives t=s+\texttt{Dt}; or
- if `Dt` is also omitted, `.default_discrimination_horizon()` uses the
  fitted maximum follow-up time and otherwise the maximum finite
  observed event or censoring time.

Every horizon must be strictly greater than its landmark. A scalar
horizon is repeated over several landmarks; otherwise the number of
horizons must equal the number of landmarks.

### 6.2 Step 2: construct a landmark-specific risk set

`.dynamic_discrimination_risk_set()` uses `.subject_event_outcomes()` to
obtain one terminal outcome per subject and `.subject_time_map()` to
align the landmark and horizon declarations.

It restricts every longitudinal history to measurements at or before the
landmark. Subjects whose outcome was already observed at or before the
landmark are subsequently removed because conditional prediction assumes
survival beyond s.

### 6.3 Step 3: predict at the exact horizon

`.discrimination_time_grid()` creates 50 equally spaced prediction times
from the landmark to the horizon, including both endpoints. Prediction
uses the same
[`predict.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/predict.JoiNMeFit.md)
event process as concordance. The row at the exact horizon is selected
using a relative floating-point tolerance; a neighbouring grid point is
not substituted.

For every subject, the risk-set data then contain:

- `risk`, equal to one minus posterior mean survival at the horizon;
- `event_time`, the observed terminal time;
- `event_status` and `event_type`;
- `event_window`, the case indicator;
- `time_start` and `time_horizon`; and
- `time_window`, observed or truncated follow-up from the landmark.

### 6.4 Step 4: classify observed cases and controls

`.auc_from_risk_set()` selects finite risks for:

- cases with `event_window == 1`; and
- controls with `event_window == 0` and `event_time > time_horizon`.

The second condition is important. It prevents a subject censored, or
experiencing a competing event, before the horizon from being
incorrectly labelled as a control.

### 6.5 Step 5: calculate and return AUC

`.auc_from_risk_set()` applies the average-rank identity above. If there
is no case or no control, no empirical case–control comparison exists
and `auc` is returned as `NA`.

The result contains:

| Column         | Meaning                                 |
|----------------|-----------------------------------------|
| `time_start`   | Landmark s                              |
| `time_horizon` | Horizon t                               |
| `auc`          | Empirical cumulative/dynamic AUC        |
| `n_cases`      | Observed cause-specific cases in (s,t\] |
| `n_controls`   | Subjects observed event-free beyond t   |
| `n_pairs`      | n\_\text{cases}n\_\text{controls}       |

## 7 Concordance and AUC answer different questions

| Feature | `concordance()` | [`auc()`](https://trinhdhk.github.io/joinme/reference/auc.md) |
|----|----|----|
| Scientific question | Are earlier failures assigned lower survival than subjects observed to survive longer? | Are cases by a specified horizon assigned greater risk than controls beyond it? |
| Time summary | One follow-up-wide index | One value per landmark–horizon pair |
| Prediction used | The survival curve at each observable event time | Risk 1-S(t\mid s) at one horizon |
| Comparison unit | Earlier-event subject versus later-observed subject | Case versus control |
| Early censoring | Pair remains usable only before censoring | Subject is omitted if status at the horizon is unknown |
| Prediction ties | Half credit | Half credit |
| Event-time ties | No strict failure order; excluded | Both may be cases and compared with controls |
| Weighting | Selectable event-time weighting | Unweighted empirical rank comparison |
| Competing event | Cause-specific censoring at its event time | Omitted if it occurs before the horizon |
| Principal output | `concordance` and weighted pair counts | `auc`, case count, control count, pair count |

Neither measure should be substituted for the other. AUC can vary
markedly with the selected horizon even when follow-up-wide concordance
changes little. Concordance can conceal clinically important
horizon-specific behaviour by averaging comparisons over event times.

## 8 Reproducible use

### 8.1 Concordance

Code

``` r

# Subject-specific origin at the final pre-outcome measurement.
concordance(fit)

# Common study-time landmark and unit pair weights.
concordance(
  fit,
  time_start = 2,
  type_weights = "none",
  n_samples = 500,
  seed = 2026
)

# Uno-style censoring event weights.
concordance(
  fit,
  time_start = 2,
  type_weights = "n/G2",
  n_samples = 500,
  seed = 2026
)
```

The default `n_samples = 200` is intended to control prediction cost.
Increasing it reduces Monte Carlo variation in posterior mean curve
ordinates. Stability can be checked by repeating the calculation with a
larger value while keeping `seed` fixed for reproducibility.

### 8.2 Time-dependent AUC

Code

``` r

# Risk at time 5 among subjects under observation after time 2.
auc(
  fit,
  time_start = 2,
  time_horizon = 5,
  n_samples = 500,
  seed = 2026
)

# A three-year prediction window from several landmarks.
auc(
  fit,
  time_start = c(0, 1, 2),
  Dt = 3,
  n_samples = 500,
  seed = 2026
)

# Equivalent public alias.
AUC(fit, time_start = 2, time_horizon = 5)
```

For several landmarks, each row is a separate cumulative/dynamic
estimand. It is not a repeated-measures sample from one common AUC.

### 8.3 Evaluation on new subjects

Code

``` r

cIndex <- concordance(
  fit,
  newdataLong = validation_long,
  newdataEvent = validation_event,
  time_start = 2,
  n_samples = 500,
  seed = 2026
)

validation_auc <- auc(
  fit,
  newdataLong = validation_long,
  newdataEvent = validation_event,
  time_start = 2,
  time_horizon = 5,
  n_samples = 500,
  seed = 2026
)
```

When new data are omitted, the fitting data are used. This is useful for
model checking but generally gives an optimistic description of
predictive performance because the same subjects informed the posterior.
Predictive claims should preferably use external validation data or a
pre-specified cross-validation procedure in which model fitting is
repeated within each training fold.

## 9 Interpretation and reporting

### 9.1 What should be reported for concordance

A manuscript should state:

1.  whether the origin was the default last pre-outcome measurement, a
    common landmark, or a subject-specific clinical landmark;
2.  the event cause and treatment of competing events;
3.  the weighting rule;
4.  whether evaluation was apparent, internally validated, or externally
    validated;
5.  `n_samples` and the random seed;
6.  the number of retained subjects and cause-specific events;
7.  the concordance estimate and weighted comparable-pair count; and
8.  that the reported value is calculated from posterior mean curves and
    has no uncertainty interval unless one was obtained separately.

An example methods statement is:

> We assessed follow-up-wide discrimination at a common landmark of two
> years using conditional survival-curve concordance. For each observed
> cause-specific event, both subjects’ posterior mean survival curves
> were compared at the earlier residual event time. Subjects censored
> before that time were not comparable, prediction ties received half
> credit, and event-time weights followed the pre-specified `n/G2` rule.

### 9.2 What should be reported for AUC

A manuscript should state:

1.  the landmark and horizon on the original time scale;
2.  the cumulative/dynamic case and control definitions;
3.  the event cause and treatment of competing events;
4.  that early-censored subjects were omitted and no IPCW was applied;
5.  the numbers of cases, controls, and case–control pairs;
6.  whether evaluation was apparent, internally validated, or externally
    validated; and
7.  `n_samples` and the random seed.

An example methods statement is:

> At a landmark of two years, we calculated the empirical
> cumulative/dynamic AUC for the five-year horizon. Cases experienced
> the cause of interest during years two to five; controls were observed
> event-free beyond year five. Subjects censored or experiencing a
> competing event before the horizon were omitted. The estimator
> compared posterior mean conditional risks, assigned half credit to
> ties, and did not use inverse-probability-of-censoring weights.

## 10 Assumptions, limitations, and diagnostic checks

### 10.1 Time scale and landmark comparability

All event times, longitudinal times, landmarks, and horizons must use
the same original time scale. A subject-specific concordance origin
defines a remaining-lifetime estimand; a common origin defines a
common-landmark estimand. These should not be combined under one
interpretation.

### 10.2 Availability of longitudinal history

A dynamic prediction requires longitudinal information observed by its
landmark. Removing subjects without such history changes the analysed
population. The returned `n_subjects` should therefore be checked
against the intended cohort.

### 10.3 Censoring

Concordance uses only observable orderings and optionally changes their
event-time weights. AUC omits subjects with unknown horizon status.
Neither procedure makes unexplained informative censoring harmless.
Sensitivity analyses should be considered when dropout is associated
with unobserved health deterioration.

### 10.4 Small comparison sets

Concordance requires at least one observable event–comparator pair. AUC
requires at least one case and one control. An `NA` estimate with zero
pairs is an undefined empirical estimand for the supplied data and
times, not a numerical failure to be replaced by zero.

### 10.5 Monte Carlo stability

The posterior mean predictions depend on the selected posterior draws.
Compare results at increasing `n_samples`, using a fixed seed, until
scientifically material digits are stable. This check concerns Monte
Carlo approximation only; it does not quantify sampling, validation, or
posterior uncertainty in the discrimination statistic.

### 10.6 Apparent versus predictive performance

Calculating either measure on the fitting data assesses apparent
discrimination. Random effects and association parameters have been
informed by those subjects, so the estimate is not an unbiased
assessment for future subjects. Held-out evaluation must preserve the
full modelling procedure, including model selection and tuning, within
the training data.

### 10.7 Discrimination is not calibration or utility

Neither statistic verifies absolute-risk calibration,
prediction-interval coverage, net benefit, or a clinically acceptable
decision threshold. Discrimination should be reported alongside
calibration and predictive uncertainty appropriate to the intended use.

## 11 Function-by-function reference

The complete discrimination path is summarised below. This table is also
a map for readers auditing the implementation in `R/diagnosis.R`.

| Function | Role |
|----|----|
| [`concordance.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/concordance.JoiNMeFit.md) | Public concordance method; validates, predicts, compares, and returns one follow-up-wide result |
| [`auc()`](https://trinhdhk.github.io/joinme/reference/auc.md) | S3 generic for time-dependent AUC |
| [`auc.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/auc.md) | Public JoiNMe AUC method; validates times and evaluates each landmark–horizon pair |
| `AUC()` | Alias of [`auc()`](https://trinhdhk.github.io/joinme/reference/auc.md) |
| `.resolve_train_data()` | Selects stored fitting data or requires a complete pair of new longitudinal and event data |
| `.constant_finite_value()` | Verifies that a subject-level time column is finite and constant within subject |
| `.discrimination_time_grid()` | Creates the 50-point landmark-to-horizon grid used by AUC prediction |
| `.default_discrimination_horizon()` | Finds a fitted or observed maximum time when no AUC horizon is supplied |
| `.concordance_time_weight()` | Validates concordance weighting and maps `"none"` to `"n"` |
| `.subject_time_map()` | Aligns scalar, vector, named-vector, or event-column time declarations by subject |
| `.subject_event_outcomes()` | Reduces right-censored or counting-process data to one terminal outcome per subject |
| `.last_preoutcome_measurement()` | Finds the default subject-specific concordance origin without using an outcome-time measurement |
| `.concordance_prediction_grid()` | Adds exact residual event times to each subject’s regular prediction grid |
| `.concordance_survival_curves()` | Restricts histories, calls prediction, and returns aligned outcomes and curves |
| `.concordance_survival_matrix()` | Aligns posterior mean curves on distinct residual event times without extrapolation |
| `.concordance_event_weights()` | Obtains event-time weights from `survival` and converts them to pair weights |
| `.concordance_from_survival_curves()` | Counts each observable concordant, discordant, or tied pair and forms \widehat C |
| `.dynamic_discrimination_risk_set()` | Builds one AUC risk set and extracts posterior mean risk at the exact horizon |
| `.auc_from_risk_set()` | Applies the Mann–Whitney rank identity and returns AUC with its effective sample counts |
| [`predict.JoiNMeFit()`](https://trinhdhk.github.io/joinme/reference/predict.JoiNMeFit.md) | Supplies conditional survival predictions from the complete fitted joint model |

## 12 References

Antolini L, Boracchi P, Biganzoli E (2005). A time-dependent
discrimination index for survival data. *Statistics in Medicine*, 24,
3927–3944. [doi:10.1002/sim.2427](https://doi.org/10.1002/sim.2427).

Rizopoulos D (2011). Dynamic predictions and prospective accuracy in
joint models for longitudinal and time-to-event data. *Biometrics*, 67,
819–829. [Publisher
record](https://onlinelibrary.wiley.com/doi/10.1111/j.1541-0420.2010.01546.x).

Therneau TM.
[`survival::concordance`](https://rdrr.io/pkg/survival/man/concordance.html):
concordance for survival and regression models. The official
documentation describes risk-set comparability, event ties, and the
available time-weighting rules. [R
documentation](https://stat.ethz.ch/R-manual/R-devel/library/survival/html/concordance.html).

Therneau TM. *Concordance* vignette for the `survival` package. The
vignette gives the computational and censoring conventions used by the
event-time weighting engine. [Official
vignette](https://stat.ethz.ch/R-manual/R-devel/library/survival/doc/concordance.pdf).
