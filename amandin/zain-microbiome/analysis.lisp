;;; Install Packages for program
(ql:quickload :lisp-stat)
(ql:quickload :data-frame)
(in-package :ls-user)

;;; Import Data - Abundance table + Metadata
(defparameter abundance-table (csv-to-data-frame #P"./data/mm_16s_hiseqs_abundance_table.csv"))



;;; Preliminary Statistics
(defun sparcity (abundance-table))
(defun shape (abundance-table))

;;; Spectral Clustering
(defun abundance-table->adjacency-table (abundance-table))
(defun adjacency-table->kernel (adjacency-table))
(defun similarity->laplacian (similarity))
