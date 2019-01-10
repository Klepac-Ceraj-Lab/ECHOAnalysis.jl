# Questions about Filemaker pro database

## Questions

- Each fecal sample has a `timepoint` associated with it in the sample ID,
  but they don't appear in the FilemakerPro database.
  How are the associated timepoints determined?
  (this is collectionNum)
- Many timepoints don't have associated dates.
  Related to the preceding question, subject 426 as an example
  has stool samples linked to timepoint 5, but in the database there's no date on tp 5.
- Are "collection#"s associated with fecal samples equivalent to timepoints?
  - YES
- Same question for lead levels "testNumber"
  - NO
- What are the `timestamp` and `timestamp_acct` fields?
  They seems to be mostly empty.
  - Associated with editing database - ignore  
- Lots of babies have missing `dueDate`s in the `Delivery` section.
  How are the gestational ages / corrected ages calculated?
  - Associated with genstatioal Age field (if no gest. age then based on birthday)

## Irregularities

These are in order of most challenging from a computational perspective to least

- StudyID 110 has a timepoint `2.5`, while all the others are integers
  - Can everything be left as integers?
    - Filter it!
- studyID 66 has two fecal sample collections listed under "collection 6" with
  different dates.
  - If these are different collections, can they have different collection numbers?
  - Since collecection nums are == timepoint, have to have same collection
  - Might be from different poops
- Under delivery method, several subjects have "multiple pregnancies" and then
  vaginal or cesarean.
  - Is it possible to stick to just "vaginal" and "cesarean"
    and perhaps record the multiple pregnancies elsewhere?
    - yes, Jen is fixing!
- Some headers have odd symbols like `#`
  - Sticking to alphanumeric is better (eg NumberOfSamples or NumSamples)
- Most parent fields are pascal case, eg `FecalSampleCollection`, but
  a few have spaces eg `Fecal with Ethanol`
  - Removing spaces is better
- Most child fields are camel case, eg `adoptionAgeMonths`,
  but some are all lower, eg `finishassessment`
  - this isn't really an issue, just odd


## Weirdness

- Fecal sample collection num == timepoint
- Mullen test num and lead test num are != timepoint (use dates)
