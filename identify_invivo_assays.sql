SELECT DISTINCT a.chembl_id as assay_chemblid, a.description as assay_description
FROM assays a 

-- First find ASSAY_organisms that are mammals by joining target_dictionary and organism_class: 
JOIN target_dictionary b ON a.assay_tax_id = b.tax_id
JOIN organism_class c ON b.tax_id = c.tax_id
-- Second find TARGET_organisms that are mammals by joining target_dictionary and organism_class: 
JOIN target_dictionary d ON a.tid = d.tid
JOIN organism_class e ON d.tax_id = e.tax_id

-- Keep assays where the BAO Ontology (BAO_0000218) is "organism-based format"
WHERE a.BAO_FORMAT = 'BAO_0000218' 
-- Keep assays where either the ASSAY_organism OR the TARGET_organism are mammals. This excludes bacteria, insects etc that are also classed as whole organisms: 
AND (c.l2 = 'Mammalia' OR  e.l2 = 'Mammalia') 
-- Exclude assay descriptions that relate to in vitro or ex vivo assays: 
AND NOT REGEXP_LIKE(lower(a.description), 'in[ -]?vitro|ex[ -]?vivo', 'i')
-- Exclude ADMET assays since these typically relate to pharmacokinetic parameters like Cmax, Tmax, Bioavailability or in vitro drug metabolism studies, and are therefore not disease or phenotypic assays:
AND a.assay_type != 'A'
-- Only include assays from published scientific literature. This excludes deposited datasets like TG-GATES that have existing good annotation. 
AND a.src_id = 1;
