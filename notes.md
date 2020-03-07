# Metadata export needs

## Basic

- `timepoint` (Int)
- `subjectID` (Int)
- `childGender` (String)
- `correctedAgeDays` (Int)
- `mother_HHS` (Int)
  - replace 9999 with missing
- `childBMI` (Float)
  - Currently calculated from `childHeight` and `childWeight` (I think from TimpointInfo, but not sure)
  - `row.childWeight / (row.childHeight^2) * 703`
- `hasBrain` (Bool)
  - white/gray/csf volume data
  - Ideally with link to source
- `hasHires` (Bool)
  - High resolution scans with segmentation done
  - Ideally, would have `hasScan` (bool) and `segmentationDone` (Bool)
    or something to know when to bug Sean/Muriel
  - Ideally with link to source
- `breastFedPercent` (Float64, 0:100)
- `birthType` (String, vaginal|cesarean)

## Medium Custom

I'm including my current julia code for these so you can see explicitly what I want.
Hopefully it's close enough to pseudo code that you can understand it.
If not, don't hesitate to ask.

- `ageLabel` from correctedAgeDays / 365
    ```julia
    if ismom(subject.sampleid)
        label="mom"
    elseif ismissing(subject.age) || isnan(subject.age)
        label=missing
    elseif subject.age <= 1
        label="1 and under"
    elseif 1 < subject.age < 2
        label="1 to 2"
    else
        label="2 and over"
    end
    ```
- `cogScore` and `cogAssessment`
  - score:
    - Composite scores for WPPSI, WISC, Mullen
    - mean of language and motor composites for Bayleys
  - assessment:
    - ["Bayleys", "Mullen", "WPPSI", "WISC"]
  - code:
    ```julia
    cogcols = [
        "languageComposite", # this is from Bayleys
        "motorComposite", # also Bayleys
        "mullen_EarlyLearningComposite",
        "fullScaleComposite", # this is from WPPSI
        "FSIQ_Composite",  # this is from WISC
        ]
    cogscores = by(allmeta, [:subject, :timepoint]) do sample
        cogs = filter(row-> in(row.metadatum, cogcols) &&
                            !ismissing(row.value),
                        sample)
        if nrow(cogs) > 0
            parent = cogs.parent_table[1]
            if parent == "Bayley's"
                nrow(cogs) != 2 && throw(ErrorException("Wrong number of rows for $parent"))
                assessment = "Bayleys"
                score = string(mean(parse.(Float64, cogs.value)))
            else
                nrow(cogs) != 1 && throw(ErrorException("Wrong number of rows for $parent"))
                metadatum = cogs.metadatum[1]
                if metadatum == "mullen_EarlyLearningComposite"
                    assessment = "Mullen"
                elseif metadatum == "fullScaleComposite"
                    assessment = "WPPSI"
                elseif metadatum == "FSIQ_Composite"
                    assessment = "WISC"
                else
                    throw(ErrorException("Couldn't figure out assessment $metadatum"))
                end
                score = cogs.value[1]
            end
        else
            assessment = missing
            score = missing
        end
        return DataFrame(metadatum    = ["cogAssessment", "cogScore"],
                         value        = [assessment,      score],
                         parent_table = ["Calculated",    "Calculated"])
    end

    dropmissing!(cogscores)
    ```
- `breastfeeding` (String)
  - categorical variable based on `breastFedPercent`
    ```julia
    map(eachrow(allmeta)) do row
        ismissing(row.breastFedPercent) && return missing
        !(row.breastFedPercent isa Number) && error(":breastFedPercent should be a number or missing")
        if row.breastFedPercent < 5
            return "exclussive formula"
        elseif row.breastFedPercent > 80
            return "exclussive breast"
        else
            return "mixed"
        end
    end

## Complicated Custom

### Medication and supplement data from preganancy (meds, abx, probx)
  
- Pegged to both mothers and kids
  - That is, all child sample rows should have this info from mothers
- Need to deal with trimesters

I'm not 100% sure, but I think the best way to deal with thisinformation would be sort of like a nested dictionary (JSON?)

So, something like

```json
pregnancyMeds: {
    zantac: {
        trimesters: [1,2],
        dose: "30mg",
        type: "acid blocker"
    },
    zoloft: {
        trimesters: [3],
        dose: ??,
        type: "SSRI"
    }
    # ...
}
pregnancyProbiotics: {...}
```

Since moms can have multiple meds
and I'd rather not have an arbitrary number of columns
that change as the data changes.

### Child meds/supplements

Same as as above

### Diet

This is a low-priority one, and I have no idea what the underlying data looks like. 
But I'd like to be in on discussions with eg Monique's gruoup about what they need if possible.