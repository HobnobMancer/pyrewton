# Compiling and Querying an Local CAZome database

## SQL Queries for finding engineering candidates

### Filter by protein mass

Filter for proteins less than 70kDa.
```sql
SELECT genbank_accession
FROM Proteins
Where Proteins.mass < 71000
```
Returns 19384 proteins

## Filter by Mass and GT removal

Filter for proteins less than 71kDa and do not contain a GT domain.
```sql
WITH GtCazymes (gt_cazymes) AS (
	SELECT DISTINCT genbank_accession
	FROM Proteins
	INNER JOIN Domains ON Proteins.protein_id = Domains.protein_id
	INNER JOIN CazyFamilies ON Domains.family_id = CazyFamilies.family_id
	Where CazyFamilies.family LIKE 'GT%'
)
SELECT genbank_accession
FROM Proteins
LEFT JOIN GtCazymes ON Proteins.genbank_accession = GtCazymes.gt_cazymes
Where (Proteins.mass < 71000) AND
	(Proteins.genbank_accession NOT IN GtCazymes)
```
Returns 15872 proteins.

## CAZy annotated

Retrieve all proteins catalogued within CAZy.
```sql
SELECT DISTINCT genbank_accession
FROM Proteins
INNER JOIN Domains ON Proteins.protein_id = Domains.protein_id
INNER JOIN CazyFamilies ON Domains.family_id = CazyFamilies.family_id
INNER JOIN Classifiers ON Domains.classifier_id = Classifiers.classifier_id
Where (Proteins.mass < 71000)AND
(CazyFamilies.family NOT LIKE 'GT%') AND 
(Classifiers.classifier == 'CAZy')
```
Returns 1284 proteins

## At least two dbCAN tools (irrespective of CAZy)

Retrieve proteins for which at least two dbCAN tools classified the protein as a CAZyme, irrespective if the protein is catalgoued in CAZy.
```sql
SELECT COUNT(DISTINCT classifier_id), Proteins.genbank_accession
FROM Domains
INNER JOIN Proteins ON Domains.protein_id = Proteins.protein_id
GROUP BY Domains.protein_id
HAVING COUNT(DISTINCT classifier_id) > 1
ORDER BY COUNT(DISTINCT classifier_id) ASC
```

## Mass, GT and at least 2 tools:

Filter:
- Mass < 71kDa
- Contains no GT domains
- At least two dbCAN tools

```sql
WITH GtCazymes (gt_cazymes) AS (
	SELECT DISTINCT genbank_accession
	FROM Proteins
	INNER JOIN Domains ON Proteins.protein_id = Domains.protein_id
	INNER JOIN CazyFamilies ON Domains.family_id = CazyFamilies.family_id
	Where CazyFamilies.family LIKE 'GT%'
), ClassifierFilter (classifier_count, class_gbk_accs) AS (
	SELECT COUNT(DISTINCT classifier_id), Proteins.genbank_accession
	FROM Domains
	INNER JOIN Proteins ON Domains.protein_id = Proteins.protein_id
	GROUP BY Domains.protein_id
	HAVING COUNT(DISTINCT classifier_id) > 1
	ORDER BY COUNT(DISTINCT classifier_id) ASC
), MassFilter (mass_gbk_accs) AS (
	SELECT DISTINCT genbank_accession
	FROM Proteins
	INNER JOIN Domains ON Proteins.protein_id = Domains.protein_id
	INNER JOIN CazyFamilies ON Domains.family_id = CazyFamilies.family_id
	INNER JOIN Classifiers ON Domains.classifier_id = Classifiers.classifier_id
	Where (Proteins.mass < 71000)
)
SELECT genbank_accession, Assemblies.assembly_accession, Taxonomies.genus, Taxonomies.species, ClassifierFilter.classifier_count
FROM Proteins
INNER JOIN Assemblies ON Proteins.assembly_id = Assemblies.assembly_id
INNER JOIN Taxonomies ON Assemblies.taxonomy_id = Taxonomies.taxonomy_id
LEFT JOIN MassFilter ON Proteins.genbank_accession = MassFilter.mass_gbk_accs
LEFT JOIN ClassifierFilter ON Proteins.genbank_accession = ClassifierFilter.class_gbk_accs
LEFT JOIN GtCazymes ON Proteins.genbank_accession = GtCazymes.gt_cazymes
WHERE (MassFilter.mass_gbk_accs is not null) AND 
	(ClassifierFilter.class_gbk_accs is not null) AND
	(Proteins.genbank_accession NOT IN GtCazymes)
ORDER BY ClassifierFilter.classifier_count ASC
```

Returns 10,318 proteins

## Mass, GT, at least 2 tools and NOT in CAZy:

Filter:
- Mass < 71kDa
- Contains no GT domains
- At least two dbCAN tools
- Not catalogued within CAZy

```sql
WITH GtCazymes (gt_cazymes) AS (
	SELECT DISTINCT genbank_accession
	FROM Proteins
	INNER JOIN Domains ON Proteins.protein_id = Domains.protein_id
	INNER JOIN CazyFamilies ON Domains.family_id = CazyFamilies.family_id
	Where CazyFamilies.family LIKE 'GT%'
), ClassifierFilter (classifier_count, class_gbk_accs) AS (
	SELECT COUNT(DISTINCT classifier_id), Proteins.genbank_accession
	FROM Domains
	INNER JOIN Proteins ON Domains.protein_id = Proteins.protein_id
	GROUP BY Domains.protein_id
	HAVING COUNT(DISTINCT classifier_id) > 1
	ORDER BY COUNT(DISTINCT classifier_id) ASC
), MassFilter (mass_gbk_accs) AS (
	SELECT DISTINCT genbank_accession
	FROM Proteins
	INNER JOIN Domains ON Proteins.protein_id = Domains.protein_id
	INNER JOIN CazyFamilies ON Domains.family_id = CazyFamilies.family_id
	INNER JOIN Classifiers ON Domains.classifier_id = Classifiers.classifier_id
	Where (Proteins.mass < 71000)
), NotCazy (not_cazy) AS (
	SELECT DISTINCT Proteins.genbank_accession
	FROM Domains
	INNER JOIN Classifiers ON Domains.classifier_id = Classifiers.classifier_id
	INNER JOIN Proteins ON Domains.protein_id = Proteins.protein_id
	WHERE Classifiers.classifier = 'CAZy'
)
SELECT genbank_accession, Assemblies.assembly_accession, Taxonomies.genus, Taxonomies.species, ClassifierFilter.classifier_count
FROM Proteins
INNER JOIN Assemblies ON Proteins.assembly_id = Assemblies.assembly_id
INNER JOIN Taxonomies ON Assemblies.taxonomy_id = Taxonomies.taxonomy_id
LEFT JOIN MassFilter ON Proteins.genbank_accession = MassFilter.mass_gbk_accs
LEFT JOIN ClassifierFilter ON Proteins.genbank_accession = ClassifierFilter.class_gbk_accs
LEFT JOIN GtCazymes ON Proteins.genbank_accession = GtCazymes.gt_cazymes
LEFT JOIN NotCazy ON Proteins.genbank_accession = NotCazy.not_cazy
WHERE (MassFilter.mass_gbk_accs is not null) AND 
	(ClassifierFilter.class_gbk_accs is not null) AND
	(Proteins.genbank_accession NOT IN GtCazymes) AND
	(Proteins.genbank_accession NOT IN NotCazy)
ORDER BY ClassifierFilter.classifier_count ASC
```

Return 9606 proteins

## Mass, GT, at least 2 tools, NOT in CAZy and has a dbCAN consensus

Filter:
- Mass < 71kDa
- Contains no GT domains
- At least two dbCAN tools
- Not catalogued within CAZy
- Has at least one dbCAN consensus CAZyme domain classification
- 
```sql
WITH GtCazymes (gt_cazymes) AS (
	SELECT DISTINCT genbank_accession
	FROM Proteins
	INNER JOIN Domains ON Proteins.protein_id = Domains.protein_id
	INNER JOIN CazyFamilies ON Domains.family_id = CazyFamilies.family_id
	Where CazyFamilies.family LIKE 'GT%'
), ClassifierFilter (classifier_count, class_gbk_accs) AS (
	SELECT COUNT(DISTINCT classifier_id), Proteins.genbank_accession
	FROM Domains
	INNER JOIN Proteins ON Domains.protein_id = Proteins.protein_id
	GROUP BY Domains.protein_id
	HAVING COUNT(DISTINCT classifier_id) > 1
	ORDER BY COUNT(DISTINCT classifier_id) ASC
), MassFilter (mass_gbk_accs) AS (
	SELECT DISTINCT genbank_accession
	FROM Proteins
	INNER JOIN Domains ON Proteins.protein_id = Domains.protein_id
	INNER JOIN CazyFamilies ON Domains.family_id = CazyFamilies.family_id
	INNER JOIN Classifiers ON Domains.classifier_id = Classifiers.classifier_id
	Where (Proteins.mass < 71000)
), NotCazy (not_cazy) AS (
	SELECT DISTINCT Proteins.genbank_accession
	FROM Domains
	INNER JOIN Classifiers ON Domains.classifier_id = Classifiers.classifier_id
	INNER JOIN Proteins ON Domains.protein_id = Proteins.protein_id
	WHERE Classifiers.classifier = 'CAZy'
), IndbCAN (in_dbcan) AS (
	SELECT DISTINCT Proteins.genbank_accession
	FROM Domains
	INNER JOIN Classifiers ON Domains.classifier_id = Classifiers.classifier_id
	INNER JOIN Proteins ON Domains.protein_id = Proteins.protein_id
	WHERE Classifiers.classifier = 'dbCAN'
)
SELECT genbank_accession, Assemblies.assembly_accession, Taxonomies.genus, Taxonomies.species, ClassifierFilter.classifier_count
FROM Proteins
INNER JOIN Assemblies ON Proteins.assembly_id = Assemblies.assembly_id
INNER JOIN Taxonomies ON Assemblies.taxonomy_id = Taxonomies.taxonomy_id
LEFT JOIN MassFilter ON Proteins.genbank_accession = MassFilter.mass_gbk_accs
LEFT JOIN ClassifierFilter ON Proteins.genbank_accession = ClassifierFilter.class_gbk_accs
LEFT JOIN GtCazymes ON Proteins.genbank_accession = GtCazymes.gt_cazymes
LEFT JOIN NotCazy ON Proteins.genbank_accession = NotCazy.not_cazy
LEFT JOIN IndbCAN ON Proteins.genbank_accession = IndbCAN.in_dbcan
WHERE (MassFilter.mass_gbk_accs is not null) AND 
	(ClassifierFilter.class_gbk_accs is not null) AND
	(Proteins.genbank_accession NOT IN GtCazymes) AND
	(Proteins.genbank_accession NOT IN NotCazy) AND 
	(Proteins.genbank_accession IN IndbCAN)
ORDER BY ClassifierFilter.classifier_count ASC
```
Return 9556 proteins

## Confused predictions

Confused predictions are defined as those with at least 2 classifiers annotating the protein as a CAZyme 
but no consensus dbCAN prediction was achieved.

Filters:
- No dbCAN consensus
- Not in CAZy
- Less than 71kDa
- Has no GT domain

```sql
WITH GtCazymes (gt_cazymes) AS (
	SELECT DISTINCT genbank_accession
	FROM Proteins
	INNER JOIN Domains ON Proteins.protein_id = Domains.protein_id
	INNER JOIN CazyFamilies ON Domains.family_id = CazyFamilies.family_id
	Where CazyFamilies.family LIKE 'GT%'
), ClassifierFilter (classifier_count, class_gbk_accs) AS (
	SELECT COUNT(DISTINCT classifier_id), Proteins.genbank_accession
	FROM Domains
	INNER JOIN Proteins ON Domains.protein_id = Proteins.protein_id
	GROUP BY Domains.protein_id
	HAVING COUNT(DISTINCT classifier_id) > 1
	ORDER BY COUNT(DISTINCT classifier_id) ASC
), MassFilter (mass_gbk_accs) AS (
	SELECT DISTINCT genbank_accession
	FROM Proteins
	INNER JOIN Domains ON Proteins.protein_id = Domains.protein_id
	INNER JOIN CazyFamilies ON Domains.family_id = CazyFamilies.family_id
	INNER JOIN Classifiers ON Domains.classifier_id = Classifiers.classifier_id
	Where (Proteins.mass < 71000)
), NotCazy (not_cazy) AS (
	SELECT DISTINCT Proteins.genbank_accession
	FROM Domains
	INNER JOIN Classifiers ON Domains.classifier_id = Classifiers.classifier_id
	INNER JOIN Proteins ON Domains.protein_id = Proteins.protein_id
	WHERE Classifiers.classifier = 'CAZy'
), IndbCAN (in_dbcan) AS (
	SELECT DISTINCT Proteins.genbank_accession
	FROM Domains
	INNER JOIN Classifiers ON Domains.classifier_id = Classifiers.classifier_id
	INNER JOIN Proteins ON Domains.protein_id = Proteins.protein_id
	WHERE Classifiers.classifier = 'dbCAN'
)
SELECT genbank_accession, Assemblies.assembly_accession, Taxonomies.genus, Taxonomies.species, ClassifierFilter.classifier_count
FROM Proteins
INNER JOIN Assemblies ON Proteins.assembly_id = Assemblies.assembly_id
INNER JOIN Taxonomies ON Assemblies.taxonomy_id = Taxonomies.taxonomy_id
LEFT JOIN MassFilter ON Proteins.genbank_accession = MassFilter.mass_gbk_accs
LEFT JOIN ClassifierFilter ON Proteins.genbank_accession = ClassifierFilter.class_gbk_accs
LEFT JOIN GtCazymes ON Proteins.genbank_accession = GtCazymes.gt_cazymes
LEFT JOIN NotCazy ON Proteins.genbank_accession = NotCazy.not_cazy
LEFT JOIN IndbCAN ON Proteins.genbank_accession = IndbCAN.in_dbcan
WHERE (MassFilter.mass_gbk_accs is not null) AND 
	(ClassifierFilter.class_gbk_accs is not null) AND
	(Proteins.genbank_accession NOT IN GtCazymes) AND
	(Proteins.genbank_accession NOT IN NotCazy) AND 
	(Proteins.genbank_accession NOT IN IndbCAN)
ORDER BY ClassifierFilter.classifier_count ASC
```
Returns 50 proteins.

Repeating the quering, after removing the 'not in CAZy' filter, returns 103 proteins.
```sql
WITH GtCazymes (gt_cazymes) AS (
	SELECT DISTINCT genbank_accession
	FROM Proteins
	INNER JOIN Domains ON Proteins.protein_id = Domains.protein_id
	INNER JOIN CazyFamilies ON Domains.family_id = CazyFamilies.family_id
	Where CazyFamilies.family LIKE 'GT%'
), ClassifierFilter (classifier_count, class_gbk_accs) AS (
	SELECT COUNT(DISTINCT classifier_id), Proteins.genbank_accession
	FROM Domains
	INNER JOIN Proteins ON Domains.protein_id = Proteins.protein_id
	GROUP BY Domains.protein_id
	HAVING COUNT(DISTINCT classifier_id) > 1
	ORDER BY COUNT(DISTINCT classifier_id) ASC
), MassFilter (mass_gbk_accs) AS (
	SELECT DISTINCT genbank_accession
	FROM Proteins
	INNER JOIN Domains ON Proteins.protein_id = Domains.protein_id
	INNER JOIN CazyFamilies ON Domains.family_id = CazyFamilies.family_id
	INNER JOIN Classifiers ON Domains.classifier_id = Classifiers.classifier_id
	Where (Proteins.mass < 71000)
), IndbCAN (in_dbcan) AS (
	SELECT DISTINCT Proteins.genbank_accession
	FROM Domains
	INNER JOIN Classifiers ON Domains.classifier_id = Classifiers.classifier_id
	INNER JOIN Proteins ON Domains.protein_id = Proteins.protein_id
	WHERE Classifiers.classifier = 'dbCAN'
)
SELECT genbank_accession, Assemblies.assembly_accession, Taxonomies.genus, Taxonomies.species, ClassifierFilter.classifier_count
FROM Proteins
INNER JOIN Assemblies ON Proteins.assembly_id = Assemblies.assembly_id
INNER JOIN Taxonomies ON Assemblies.taxonomy_id = Taxonomies.taxonomy_id
LEFT JOIN MassFilter ON Proteins.genbank_accession = MassFilter.mass_gbk_accs
LEFT JOIN ClassifierFilter ON Proteins.genbank_accession = ClassifierFilter.class_gbk_accs
LEFT JOIN GtCazymes ON Proteins.genbank_accession = GtCazymes.gt_cazymes
LEFT JOIN IndbCAN ON Proteins.genbank_accession = IndbCAN.in_dbcan
WHERE (MassFilter.mass_gbk_accs is not null) AND 
	(ClassifierFilter.class_gbk_accs is not null) AND
	(Proteins.genbank_accession NOT IN GtCazymes) AND 
	(Proteins.genbank_accession NOT IN IndbCAN)
ORDER BY ClassifierFilter.classifier_count ASC
```
