# Bytecode interpreter

This directory describes the small transformation language shared by the R
compiler and the R and Stan interpreters. The interpreter is deliberately
independent of joint models: it receives an input value, an integer instruction
stream, a numeric constant pool and optional instruction-aligned affine
indices. It neither reads a fitted model nor knows what an association or link
function represents.

## Programme representation

A programme is reverse Polish notation. `PUSH_X` places the supplied scalar
input on the stack. `PUSH_CONST` consumes the next constant from the constant
pool. Unary instructions replace the stack top, and binary instructions pop
the right operand before combining it with the left operand. A valid programme
finishes with exactly one stack value.

|   Code | Instruction                        | Arity | Result                                |
| -----: | ---------------------------------- | ----: | ------------------------------------- |
|      0 | `PUSH_X`                         |     0 | input value                           |
|      1 | `PUSH_CONST`                     |     0 | next constant                         |
|      2 | `ADD`                            |     2 | left + right                          |
|      3 | `SUB`                            |     2 | left - right                          |
|      4 | `MUL`                            |     2 | left × right                         |
|      5 | `DIV`                            |     2 | left ÷ right                         |
|      6 | `LOG`                            |     1 | natural logarithm                     |
|      7 | `EXP`                            |     1 | exponential                           |
|      8 | `SQRT`                           |     1 | square root                           |
|      9 | `INV_LOGIT`                      |     1 | logistic inverse link                 |
|     10 | `LOGIT`                          |     1 | logistic link                         |
|     11 | `RECIPROCAL`                     |     1 | reciprocal                            |
|     12 | `POW`                            |     2 | left raised to right                  |
| 13–16 | `SIN`, `COS`, `TAN`, `ABS` |     1 | named operation                       |
|     17 | `SQUARE`                         |     1 | square                                |
| 18–23 | hyperbolic functions and inverses  |     1 | named operation                       |
|     24 | `SOFTPLUS`                       |     1 | log(1 + exp(value))                   |
|     25 | `CBRT`                           |     1 | real cube root                        |
|     26 | `PHI`                            |     1 | standard Normal distribution function |
|     27 | `INV_PHI`                        |     1 | standard Normal quantile              |

For backward compatibility, a programme beginning with a unary instruction is
read as though `PUSH_X` preceded it.
