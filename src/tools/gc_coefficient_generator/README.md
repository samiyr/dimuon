# SIDIS coefficient generator

The Mathematica code in `sidis_coefficients.wl` includes everything needed to generate the `gc` coefficients from a FORM file, as well as to make the reference tables for testing. The FORM file `ancillary_EW.frm` is taken from [[JHEP 10 (2025) 016]](https://doi.org/10.1007/JHEP10(2025)016) and modified slightly (comments removed as the parser doesn't support that).

## Generating coefficients

As an example, the following code generates the LO coefficients:
```mma
FORMfile = "ancillary_EW.frm";
FORMcoefficients = Import[FORMfile, "Text"];
FORMinput = FORMreplacement[FORMcoefficients];
ToExpression[FORMinput]

lo = ({
    	eF1C0q2qM, eF2C0q2qM, eF3C0q2qM
    } = {
     	evaluate[f1[2 FTC0q2qM]], evaluate[f2[2 FTC0q2qM, 0]], 
    	evaluate[2 F3C0q2qM]
    }
);

loOptimized = optimize[lo, 4];

loNames = {"F1C0q2qM", "F2C0q2qM", "F3C0q2qM"};
generateCode[
	"gc/" <> # <> ".cpp" & /@ loNames,
	"gc/lo_cache.cpp", "gc/lo_decl.cpp", "gc/lo_vars.cpp",
	"lo",
	loIR
] // export;
```

Note the normalization factors of 2. The `4` in `optimize[lo, 4]` is the minimum leaf count for the common subexpression elimination; see section "`gc` coefficients" in `numerical_evaluation.md` for more. The various files that are created and exported are also explained in `numerical_evaluation.md`. The other orders work analogously. The order and naming conventions are as follows for NLO

```mma
nlo = ({
    	eF1C1q2qM, eF2C1q2qM, eF3C1q2qM,
    	eF1C1q2gM, eF2C1q2gM, eF3C1q2gM,
    	eF1C1g2qM, eF2C1g2qM, eF3C1g2qM
    } = {
    	evaluate[f1[2 FTC1q2qM]], evaluate[f2[2 FTC1q2qM, 2 FLC1q2qM]], 
     	evaluate[2 F3C1q2qM],
     	evaluate[f1[2 FTC1q2gM]], evaluate[f2[2 FTC1q2gM, 2 FLC1q2gM]], 
     	evaluate[2 F3C1q2gM],
     	evaluate[f1[2 FTC1g2qM]], evaluate[f2[2 FTC1g2qM, 2 FLC1g2qM]], 
     	evaluate[2 F3C1g2qM]
    }
);

nloNames = {
   "F1C1q2qM", "F2C1q2qM", "F3C1q2qM",
   "F1C1q2gM", "F2C1q2gM", "F3C1q2gM",
   "F1C1g2qM", "F2C1g2qM", "F3C1g2qM"
};
```

and for NNLO

```mma
nnlo = ({
     	eF1C2q2qM, eF2C2q2qM, eF3C2q2qM,
     	eF1C2q2qMFcon2, eF2C2q2qMFcon2,
     	eF1C2q2qMFcon1, eF2C2q2qMFcon1, eF3C2q2qMFcon1,
     	eF1C2q2qbM, eF2C2q2qbM, eF3C2q2qbM,
     	eF1C2q2qbMFcon, eF2C2q2qbMFcon, eF3C2q2qbMFcon,
     	eF1C2q2qpM1, eF2C2q2qpM1, eF3C2q2qpM1,
     	eF1C2q2qpM2, eF2C2q2qpM2, eF3C2q2qpM2,
     
     	eF1C2q2gM, eF2C2q2gM, eF3C2q2gM,
     	eF1C2g2qM, eF2C2g2qM, eF3C2g2qM,
     
     	eF1C2g2gM, eF2C2g2gM
    } = {
      	evaluate[f1[2 FTC2q2qM]], evaluate[f2[2 FTC2q2qM, 2 FLC2q2qM]], 
      	evaluate[2 F3C2q2qM],
      	evaluate[f1[2 FTC2q2qMFcon2]], 
      	evaluate[f2[2 FTC2q2qMFcon2, 2 FLC2q2qMFcon2]],
      	evaluate[f1[2 FTC2q2qMFcon1]], 
      	evaluate[f2[2 FTC2q2qMFcon1, 2 FLC2q2qMFcon1]], 
      	evaluate[2 F3C2q2qMFcon1],
      	evaluate[f1[2 FTC2q2qbM]], evaluate[f2[2 FTC2q2qbM, 2 FLC2q2qbM]],
       	evaluate[2 F3C2q2qbM],
      	evaluate[f1[2 FTC2q2qbMFcon]], 
      	evaluate[f2[2 FTC2q2qbMFcon, 2 FLC2q2qbMFcon]], 
      	evaluate[2 F3C2q2qbMFcon],
      	evaluate[f1[2 FTC2q2qpM1]], 
      	evaluate[f2[2 FTC2q2qpM1, 2 FLC2q2qpM1]], evaluate[2 F3C2q2qpM1],
      	evaluate[f1[2 FTC2q2qpM2]], 
      	evaluate[f2[2 FTC2q2qpM2, 2 FLC2q2qpM2]], evaluate[2 F3C2q2qpM2],
      
      	evaluate[f1[2 FTC2q2gM]], evaluate[f2[2 FTC2q2gM, 2 FLC2q2gM]], 
      	evaluate[2 F3C2q2gM],
      	evaluate[f1[2 FTC2g2qM]], evaluate[f2[2 FTC2g2qM, 2 FLC2g2qM]], 
      	evaluate[2 F3C2g2qM],
      
      	evaluate[f1[2 FTC2g2gM]], evaluate[f2[2 FTC2g2gM, 2 FLC2g2gM]]
    } // parallelize
);

nnloNames = {
	"F1C2q2qM", "F2C2q2qM", "F3C2q2qM", "F1C2q2qMFcon2", 
   	"F2C2q2qMFcon2", "F1C2q2qMFcon1", "F2C2q2qMFcon1", "F3C2q2qMFcon1",
    "F1C2q2qbM", "F2C2q2qbM", "F3C2q2qbM", "F1C2q2qbMFcon", 
   	"F2C2q2qbMFcon", "F3C2q2qbMFcon", "F1C2q2qpM1", "F2C2q2qpM1", 
   	"F3C2q2qpM1", "F1C2q2qpM2", "F2C2q2qpM2", "F3C2q2qpM2", "F1C2q2gM",
    "F2C2q2gM", "F3C2q2gM", "F1C2g2qM", "F2C2g2qM", "F3C2g2qM", 
   	"F1C2g2gM", "F2C2g2gM"
};
```

## Generating reference tables

Using the generated coefficients (`lo`, `nlo`, and `nnlo` from the previous section), the LO and NLO integrals can be computed analytically with e.g.
```mma
loIntegrals = calculateIntegralWith[#, pdfF, ffF] & /@ lo;
```
In the above, `pdfF` and `ffF` are the toy PDF and FF functions. For NLO, the command is essentially the same, although it is better to wrap it in `ParallelMap` to speed things up:
```mma
nloIntegrals = ParallelMap[calculateIntegralWith[#, pdfF, ffF] &, nlo];
```

At NNLO, however, analytical integration is no longer possible (not with Mathematica at least). Instead, they have to be computed numerically with
```mma
nnloIntegrals = setupIntegrals /@ nnlo;
```

This doesn't actually compute any integrals, just sets them up, since numerical integration cannot be parametric. The actual computation is done with
```mma
calculateIntegral[coeff, {xx, zz}, {pdfF, ffF}, {
	NC -> 3, NF -> 4, LMUF -> Log[4], LMUA -> Log[4], LMUR -> Log[4]
}, 30]
```

where `xx` and `zz` are variables defining $x$ and $z$. The replacement rules specify other parameter values. The `30` argument is the precision: machine-precision is not enough for all coefficients. Using higher precision unfortunately means that the code will take multiple hours to run.