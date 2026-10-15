package main

// HTTP handlers of the lens-group endpoints (POST /api/lens/group/*) and of
// POST /api/lens/parax. They decode the request, call the shared run* core of
// group_core.go and answer with the numeric response structures; the CLI
// subcommands of group_cmd.go call the very same functions.

import (
	"net/http"
)

// ---------------------------------------------------------------------------
// The endpoints
// ---------------------------------------------------------------------------

// handleLensGroupLayout answers POST /api/lens/group/layout: the global vertex,
// axis, normal and fold count of every surface (same core as
// `oz lens group layout`).
func handleLensGroupLayout(w http.ResponseWriter, r *http.Request) {
	var req groupRequest
	if !decodePOST(w, r, &req, 4<<20) {
		return
	}
	res, err := runGroupLayout(&req)
	if err != nil {
		writeAPIErr(w, http.StatusBadRequest, err)
		return
	}
	writeJSON(w, http.StatusOK, res)
}

// handleLensGroupImage answers POST /api/lens/group/image: the best-focus image
// of the bundle through the group (same core as `oz lens group image`).
func handleLensGroupImage(w http.ResponseWriter, r *http.Request) {
	var req groupRequest
	if !decodePOST(w, r, &req, 4<<20) {
		return
	}
	res, err := runGroupImage(&req)
	if err != nil {
		writeAPIErr(w, http.StatusBadRequest, err)
		return
	}
	writeJSON(w, http.StatusOK, res)
}

// handleLensGroupAfocal answers POST /api/lens/group/afocal: the angular
// magnification Γ of an afocal group (same core as `oz lens group afocal`).
func handleLensGroupAfocal(w http.ResponseWriter, r *http.Request) {
	var req groupRequest
	if !decodePOST(w, r, &req, 4<<20) {
		return
	}
	res, err := runGroupAfocal(&req)
	if err != nil {
		writeAPIErr(w, http.StatusBadRequest, err)
		return
	}
	writeJSON(w, http.StatusOK, res)
}

// handleLensGroupDispersion answers POST /api/lens/group/dispersion: the axial
// focus shift between two wavelengths (same core as `oz lens group dispersion`).
func handleLensGroupDispersion(w http.ResponseWriter, r *http.Request) {
	var req groupRequest
	if !decodePOST(w, r, &req, 4<<20) {
		return
	}
	res, err := runGroupDispersion(&req)
	if err != nil {
		writeAPIErr(w, http.StatusBadRequest, err)
		return
	}
	writeJSON(w, http.StatusOK, res)
}

// handleLensGroupSet answers POST /api/lens/group/set: the parameter edits are
// applied to the group in order and the updated group is returned (same core as
// `oz lens group set`; the error of an edit is the error of the API).
func handleLensGroupSet(w http.ResponseWriter, r *http.Request) {
	var req groupSetRequest
	if !decodePOST(w, r, &req, 4<<20) {
		return
	}
	res, err := runGroupSet(&req)
	if err != nil {
		writeAPIErr(w, http.StatusBadRequest, err)
		return
	}
	writeJSON(w, http.StatusOK, res)
}

// handleLensGroupPrescription answers POST /api/lens/group/prescription: the
// card form of the group (same core as `oz lens group to-prescription`).
func handleLensGroupPrescription(w http.ResponseWriter, r *http.Request) {
	var req groupRequest
	if !decodePOST(w, r, &req, 4<<20) {
		return
	}
	p, err := runGroupPrescription(&req)
	if err != nil {
		writeAPIErr(w, http.StatusBadRequest, err)
		return
	}
	writeJSON(w, http.StatusOK, map[string]any{"prescription": p})
}

// handleLensParax answers POST /api/lens/parax: the first-order quantities of a
// prescription, mirrors included (same core as `oz lens parax`).
func handleLensParax(w http.ResponseWriter, r *http.Request) {
	var req paraxRequest
	if !decodePOST(w, r, &req, 4<<20) {
		return
	}
	res, err := runParax(&req)
	if err != nil {
		writeAPIErr(w, http.StatusBadRequest, err)
		return
	}
	writeJSON(w, http.StatusOK, res)
}
