/*
 * Generated by asn1c-0.9.29 (http://lionet.info/asn1c)
 * From ASN.1 module "SUPL-POS"
 * 	found in "../ulp.asn1"
 * 	`asn1c -S ../../skeletons -pdu=ULP-PDU -pdu=SUPLINIT -fcompound-names -no-gen-OER`
 */

#include "SUPLPOS.h"

asn_TYPE_member_t asn_MBR_SUPLPOS_1[] = {
	{ ATF_NOFLAGS, 0, offsetof(struct SUPLPOS, posPayLoad),
		(ASN_TAG_CLASS_CONTEXT | (0 << 2)),
		+1,	/* EXPLICIT tag at current level */
		&asn_DEF_PosPayLoad,
		0,
		{ 0, 0, 0 },
		0, 0, /* No default value */
		"posPayLoad"
		},
	{ ATF_POINTER, 2, offsetof(struct SUPLPOS, velocity),
		(ASN_TAG_CLASS_CONTEXT | (1 << 2)),
		+1,	/* EXPLICIT tag at current level */
		&asn_DEF_Velocity,
		0,
		{ 0, 0, 0 },
		0, 0, /* No default value */
		"velocity"
		},
	{ ATF_POINTER, 1, offsetof(struct SUPLPOS, ver2_SUPL_POS_extension),
		(ASN_TAG_CLASS_CONTEXT | (2 << 2)),
		-1,	/* IMPLICIT tag at current level */
		&asn_DEF_Ver2_SUPL_POS_extension,
		0,
		{ 0, 0, 0 },
		0, 0, /* No default value */
		"ver2-SUPL-POS-extension"
		},
};
static const int asn_MAP_SUPLPOS_oms_1[] = { 1, 2 };
static const ber_tlv_tag_t asn_DEF_SUPLPOS_tags_1[] = {
	(ASN_TAG_CLASS_UNIVERSAL | (16 << 2))
};
static const asn_TYPE_tag2member_t asn_MAP_SUPLPOS_tag2el_1[] = {
    { (ASN_TAG_CLASS_CONTEXT | (0 << 2)), 0, 0, 0 }, /* posPayLoad */
    { (ASN_TAG_CLASS_CONTEXT | (1 << 2)), 1, 0, 0 }, /* velocity */
    { (ASN_TAG_CLASS_CONTEXT | (2 << 2)), 2, 0, 0 } /* ver2-SUPL-POS-extension */
};
asn_SEQUENCE_specifics_t asn_SPC_SUPLPOS_specs_1 = {
	sizeof(struct SUPLPOS),
	offsetof(struct SUPLPOS, _asn_ctx),
	asn_MAP_SUPLPOS_tag2el_1,
	3,	/* Count of tags in the map */
	asn_MAP_SUPLPOS_oms_1,	/* Optional members */
	1, 1,	/* Root/Additions */
	2,	/* First extension addition */
};
asn_TYPE_descriptor_t asn_DEF_SUPLPOS = {
	"SUPLPOS",
	"SUPLPOS",
	&asn_OP_SEQUENCE,
	asn_DEF_SUPLPOS_tags_1,
	sizeof(asn_DEF_SUPLPOS_tags_1)
		/sizeof(asn_DEF_SUPLPOS_tags_1[0]), /* 1 */
	asn_DEF_SUPLPOS_tags_1,	/* Same as above */
	sizeof(asn_DEF_SUPLPOS_tags_1)
		/sizeof(asn_DEF_SUPLPOS_tags_1[0]), /* 1 */
	{ 0, 0, SEQUENCE_constraint },
	asn_MBR_SUPLPOS_1,
	3,	/* Elements count */
	&asn_SPC_SUPLPOS_specs_1	/* Additional specs */
};
