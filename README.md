# Physics modeling in Rust

This repo contains code for simulating electrons flying around a biconic cusp.

Terms:
- electrons are the tiny charged particles we associate with atoms
- a cusp is a section of narrowing magnetic fields
- a biconic cusp is a regionw where two magnetic fields tighten in opposite directions

The goal is to provide a computation framework for modeling nuclear fusion plasmas using particle-in-cell techniques. So far, I've implemented the "particle" but we have not yet implemented the "cell".

## Purpose

The primary purpose of this work is to evaluate the usage of Rust for plasma physics modeling and simulation.

Typically, one would reach for Python, as [I did previously](https://github.com/nuclearfusion2017/electron-optimization/blob/master/potential_optimizer.py), but you will eventually need to turn to a lower level language for the real performance sensitive stuff.

## Current evaluation

- Scipy is missing, though I have started on building parts of it at http://github.com/jkelleyrtp/integrals
- Still no great plotting tool
- Very performant!
- No way to parallelize onto GPU using native rust code
- Math is a tricky thing to translate
- Builtin testing framework is very nice
- Datatype conversions being explicit is kinda okay
- An offset system would be nice too for better accuracy with lower precision
-
