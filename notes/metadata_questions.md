# Questions about Filemaker pro database

## Questions

- Each fecal sample has a `timepoint` associated with it in the sample ID,
  but they don't appear in the FilemakerPro database.
  How are the associated timepoints determined?
- What are the `timestamp` and `timestamp_acct` fields?
  They seems to be mostly empty.
- Lots of babies have missing `dueDate`s in the `Delivery` section.
  How are the gestational ages / corrected ages calculated?


## Irregularities

These are in order of most challenging from a computational perspective to least

- StudyID 110 has a timepoint `2.5`, while all the others are integers
  - Can everything be left as integers?
- studyID 66 has two fecal sample collections listed under "collection 6" with
  different dates.
  - If these are different collections, can they have different collection numbers?
- Some headers have odd symbols like `#`
  - Sticking to alphanumeric is better (eg NumberOfSamples or NumSamples)
- Most parent fields are pascal case, eg `FecalSampleCollection`, but
  a few have spaces eg `Fecal with Ethanol`
  - Removing spaces is better
- Most child fields are camel case, eg `adoptionAgeMonths`,
  but some are all lower, eg `finishassessment`
  - this isn't really an issue, just odd
