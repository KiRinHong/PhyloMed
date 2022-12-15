#!/bin/bash
rm medtest.dag.*
rm -r GETSIM/logs/
mkdir GETSIM/logs/
rm COMBSIM/logs/*
rm OUTPUT/logs/*
rm -r COMBSIM/shared/Type*/
mkdir COMBSIM/shared/Type{1..2}

