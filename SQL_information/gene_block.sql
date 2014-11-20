-- MySQL dump 10.13  Distrib 5.5.37, for debian-linux-gnu (x86_64)
--
-- Host: localhost    Database: gene_block
-- ------------------------------------------------------
-- Server version	5.5.37-0ubuntu0.13.10.1

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `blast_hit`
--

DROP TABLE IF EXISTS `blast_hit`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `blast_hit` (
  `hit_id` int(11) NOT NULL AUTO_INCREMENT,
  `query_locus` varchar(45) NOT NULL,
  `subject_locus` varchar(45) NOT NULL,
  `aligned_length` int(11) NOT NULL,
  `bits_score` float NOT NULL,
  `e_val` float NOT NULL,
  `align_subject_stop` int(11) NOT NULL,
  `align_subject_start` int(11) NOT NULL,
  `align_query_stop` int(11) NOT NULL,
  `align_query_start` int(11) NOT NULL,
  `number_gaps` int(11) NOT NULL,
  `number_mismatched` int(11) NOT NULL,
  `percent_ident` float NOT NULL,
  PRIMARY KEY (`hit_id`),
  UNIQUE KEY `hit_id_UNIQUE` (`hit_id`),
  KEY `locus_idx` (`subject_locus`),
  KEY `blast_query_locus_idx` (`query_locus`),
  CONSTRAINT `locus` FOREIGN KEY (`subject_locus`) REFERENCES `orf_entry` (`locus`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `blast_query_locus` FOREIGN KEY (`query_locus`) REFERENCES `blast_query` (`blast_query_locus`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COMMENT='This table sotres the result of a blast hit.  It contains information about the  BLAST hit, as well as the target and query sequences.';
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `blast_hit`
--

LOCK TABLES `blast_hit` WRITE;
/*!40000 ALTER TABLE `blast_hit` DISABLE KEYS */;
/*!40000 ALTER TABLE `blast_hit` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `blast_query`
--

DROP TABLE IF EXISTS `blast_query`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `blast_query` (
  `blast_query_locus` varchar(45) NOT NULL,
  `accession_no` varchar(45) NOT NULL,
  `start` int(11) NOT NULL,
  `stop` int(11) NOT NULL,
  `strand` int(11) NOT NULL COMMENT 'may only contain 1, -1\n',
  `genbank_annotation` varchar(10) DEFAULT NULL,
  `product_type` varchar(45) DEFAULT NULL,
  `translation_table` varchar(45) DEFAULT NULL,
  `dna_sequence` mediumtext,
  `amino_acid_sequence` mediumtext,
  `protein_id` varchar(45) DEFAULT NULL,
  PRIMARY KEY (`blast_query_locus`),
  KEY `accession_no_idx` (`accession_no`),
  CONSTRAINT `accession_no3` FOREIGN KEY (`accession_no`) REFERENCES `sequence_info` (`accession_no`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `locus1` FOREIGN KEY (`blast_query_locus`) REFERENCES `orf_entry` (`locus`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COMMENT='This is a table that shares a schema with orf_entry. The only difference is that this table contains only sequences used in a BLAST query.';
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `blast_query`
--

LOCK TABLES `blast_query` WRITE;
/*!40000 ALTER TABLE `blast_query` DISABLE KEYS */;
/*!40000 ALTER TABLE `blast_query` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `gene_block_info`
--

DROP TABLE IF EXISTS `gene_block_info`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `gene_block_info` (
  `gene_block_name` varchar(100) NOT NULL,
  `category` varchar(45) DEFAULT NULL,
  `function` varchar(45) DEFAULT NULL,
  `notes` varchar(250) DEFAULT NULL,
  PRIMARY KEY (`gene_block_name`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COMMENT='This table will store information about a gene block.  I am currently unsure what the totality of this information is.';
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `gene_block_info`
--

LOCK TABLES `gene_block_info` WRITE;
/*!40000 ALTER TABLE `gene_block_info` DISABLE KEYS */;
/*!40000 ALTER TABLE `gene_block_info` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `gene_block_membership`
--

DROP TABLE IF EXISTS `gene_block_membership`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `gene_block_membership` (
  `gene` varchar(10) NOT NULL,
  `gene_block_name` varchar(100) DEFAULT NULL,
  PRIMARY KEY (`gene`),
  KEY `gene_block_name_idx` (`gene_block_name`),
  CONSTRAINT `gene_block_name` FOREIGN KEY (`gene_block_name`) REFERENCES `gene_block_info` (`gene_block_name`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COMMENT='This table will store information about gene_name to gene block membership. It will answer the question: given a gene annotation, which gene block does this gene belong.';
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `gene_block_membership`
--

LOCK TABLES `gene_block_membership` WRITE;
/*!40000 ALTER TABLE `gene_block_membership` DISABLE KEYS */;
/*!40000 ALTER TABLE `gene_block_membership` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `hgt`
--

DROP TABLE IF EXISTS `hgt`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `hgt` (
  `hgt_id` int(11) NOT NULL AUTO_INCREMENT,
  `accession_no` varchar(45) NOT NULL,
  `method` varchar(45) DEFAULT NULL,
  `start` int(11) DEFAULT NULL,
  `stop` int(11) DEFAULT NULL,
  `score` float DEFAULT NULL,
  PRIMARY KEY (`hgt_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COMMENT='This table contains information about HGT regions within a genomic segment, as defined by the accession number.  \nIt will be up to the query parameters to determine if a particular gene is contained within a predicted region.  More than one HGT perdiction algorithm may overlap a region.';
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `hgt`
--

LOCK TABLES `hgt` WRITE;
/*!40000 ALTER TABLE `hgt` DISABLE KEYS */;
/*!40000 ALTER TABLE `hgt` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `jobs`
--

DROP TABLE IF EXISTS `jobs`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `jobs` (
  `job_id` int(11) NOT NULL AUTO_INCREMENT,
  `user_id` int(11) DEFAULT NULL,
  `job_status` varchar(45) DEFAULT NULL,
  `active` tinyint(1) DEFAULT NULL,
  PRIMARY KEY (`job_id`),
  KEY `user_id_idx` (`user_id`),
  CONSTRAINT `user_id` FOREIGN KEY (`user_id`) REFERENCES `users` (`user_id`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COMMENT='This table will be used in job querying for the eventual web site that we are creating as part of my project.  It does not contain information about gene blocks.';
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `jobs`
--

LOCK TABLES `jobs` WRITE;
/*!40000 ALTER TABLE `jobs` DISABLE KEYS */;
/*!40000 ALTER TABLE `jobs` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `orf_entry`
--

DROP TABLE IF EXISTS `orf_entry`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `orf_entry` (
  `locus` varchar(45) NOT NULL,
  `accession_no` varchar(45) NOT NULL,
  `start` int(11) NOT NULL,
  `stop` int(11) NOT NULL,
  `strand` int(11) NOT NULL COMMENT 'may only contain 1, -1\n',
  `genbank_annotation` varchar(10) DEFAULT NULL,
  `product_type` varchar(45) DEFAULT NULL,
  `translation_table` varchar(45) DEFAULT NULL,
  `dna_sequence` mediumtext,
  `amino_acid_sequence` mediumtext,
  `protein_id` varchar(45) DEFAULT NULL,
  PRIMARY KEY (`locus`),
  KEY `accession_no_idx` (`accession_no`),
  KEY `gene_idx` (`genbank_annotation`),
  CONSTRAINT `accession_no` FOREIGN KEY (`accession_no`) REFERENCES `sequence_info` (`accession_no`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `gene` FOREIGN KEY (`genbank_annotation`) REFERENCES `gene_block_membership` (`gene`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COMMENT='This table contains detailed information about operon reading frames as annotated in genback for an organism. ';
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `orf_entry`
--

LOCK TABLES `orf_entry` WRITE;
/*!40000 ALTER TABLE `orf_entry` DISABLE KEYS */;
/*!40000 ALTER TABLE `orf_entry` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `sequence_info`
--

DROP TABLE IF EXISTS `sequence_info`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `sequence_info` (
  `accession_no` varchar(45) NOT NULL,
  `internal_tax_id` int(11) NOT NULL,
  `sequence_type` varchar(45) DEFAULT NULL,
  `version` varchar(45) DEFAULT NULL,
  `date_updated` date DEFAULT NULL,
  `sequence` longtext,
  `common_name` varchar(200) DEFAULT NULL,
  PRIMARY KEY (`accession_no`),
  KEY `tax_id_idx` (`internal_tax_id`),
  CONSTRAINT `tax_id` FOREIGN KEY (`internal_tax_id`) REFERENCES `taxonomy` (`tax_id`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `sequence_info`
--

LOCK TABLES `sequence_info` WRITE;
/*!40000 ALTER TABLE `sequence_info` DISABLE KEYS */;
/*!40000 ALTER TABLE `sequence_info` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `taxonomy`
--

DROP TABLE IF EXISTS `taxonomy`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `taxonomy` (
  `tax_id` int(11) NOT NULL AUTO_INCREMENT,
  `kingdom` varchar(45) DEFAULT NULL,
  `phylum` varchar(45) DEFAULT NULL,
  `class` varchar(45) DEFAULT NULL,
  `order` varchar(45) DEFAULT NULL,
  `family` varchar(45) DEFAULT NULL,
  `genus` varchar(45) DEFAULT NULL,
  `species` varchar(45) DEFAULT NULL,
  PRIMARY KEY (`tax_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `taxonomy`
--

LOCK TABLES `taxonomy` WRITE;
/*!40000 ALTER TABLE `taxonomy` DISABLE KEYS */;
/*!40000 ALTER TABLE `taxonomy` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `users`
--

DROP TABLE IF EXISTS `users`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `users` (
  `user_id` int(11) NOT NULL AUTO_INCREMENT,
  `user_name` varchar(45) DEFAULT NULL,
  `pwd` varchar(16) DEFAULT NULL,
  `email` varchar(62) DEFAULT NULL,
  `permission_level` varchar(45) DEFAULT NULL,
  PRIMARY KEY (`user_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COMMENT='In the case that i have users of the site, this table will contain info/ permissions.  This will need some work.';
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `users`
--

LOCK TABLES `users` WRITE;
/*!40000 ALTER TABLE `users` DISABLE KEYS */;
/*!40000 ALTER TABLE `users` ENABLE KEYS */;
UNLOCK TABLES;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2014-07-22 13:55:41
