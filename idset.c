/* PANDAseq -- Assemble paired FASTQ Illumina reads and strip the region between amplification primers.
     Copyright (C) 2011-2012  Andre Masella

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */
#include "config.h"
#include "pandaseq.h"
#include <stdlib.h>
#ifdef HAVE_PTHREAD
#        include <pthread.h>

/* All idsets share a single mutex to control reference counts */
static pthread_mutex_t ref_lock = PTHREAD_MUTEX_INITIALIZER;
#endif

struct node {
	panda_seq_identifier id;
	struct node *left;
	struct node *right;
};

struct panda_idset {
	volatile size_t refcnt;
	struct node *root;

};

PandaSet *
panda_idset_new(
	void) {
	PandaSet *set = malloc(sizeof(struct panda_idset));
	set->refcnt = 1;
	set->root = NULL;

	return set;
}

PandaSet *
panda_idset_ref(
	PandaSet * set) {
#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&ref_lock);
#endif
	set->refcnt++;
#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&ref_lock);
#endif
	return set;
}

void
node_free(
	struct node *node) {
	if (node != NULL) {
		node_free(node->left);
		node_free(node->right);
		free(node);
	}
}

void
panda_idset_unref(
	PandaSet * set) {
	int count;
#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&ref_lock);
#endif
	count = --(set->refcnt);
#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&ref_lock);
#endif
	if (count == 0) {
		node_free(set->root);
		free(set);
	}
}

void
insert(
	struct node **node,
	const panda_seq_identifier *id) {
	if (*node == NULL) {
		(*node) = malloc(sizeof(struct node));
		(*node)->left = NULL;
		(*node)->right = NULL;
		panda_seqid_clear(&(*node)->id);
		(*node)->id = *id;
	} else {
		int comparison = panda_seqid_compare(id, &(*node)->id);
		if (comparison < 0) {
			insert(&(*node)->left, id);
			if ((*node)->left->id.x < (*node)->id.x) {
				struct node *temp = (*node)->left;
				(*node)->left = temp->right;
				temp->right = *node;
				*node = temp;
			}
		} else if (comparison > 0) {
			insert(&(*node)->right, id);
			if ((*node)->right->id.x < (*node)->id.x) {
				struct node *temp = (*node)->right;
				(*node)->right = temp->left;
				temp->left = *node;
				*node = temp;
			}
		}
	}
}

void
panda_idset_add(
	PandaSet * set,
	const panda_seq_identifier *id) {
	insert(&set->root, id);
}

bool
panda_idset_add_str(
	PandaSet * set,
	const char *id,
	PandaTagging policy,
	bool *old,
	const char **end_ptr) {
	panda_seq_identifier seq_id;
	if (panda_seqid_parse_fail(&seq_id, id, policy, old, end_ptr) == 0) {
		return false;
	} else {
		panda_idset_add(set, &seq_id);
		return true;
	}
}

bool
panda_idset_contains(
	PandaSet * set,
	const panda_seq_identifier *id) {
	struct node *curr = set->root;
	while (curr != NULL) {
		int comp = panda_seqid_compare(id, &curr->id);
		if (comp == 0)
			return true;
		if (comp < 0)
			curr = curr->left;
		else
			curr = curr->right;
	}
	return false;
}
