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

#ifndef _PANDASEQ_PLUGIN_H
#        define _PANDASEQ_PLUGIN_H
#        ifdef __cplusplus
#                define EXTERN_C_BEGIN  extern "C" {
#                define EXTERN_C_END    }
#        else
#                define EXTERN_C_BEGIN
#                define EXTERN_C_END
#        endif
#        include <pandaseq.h>
#        include <stdlib.h>
#        include <string.h>
EXTERN_C_BEGIN
#        define PANDA_STRUCT_DUP(item) memcpy(malloc(sizeof(*item)), item, sizeof(*item))
#        define OPEN bool PANDACONCAT(PANDASEQ_MODULE,_LTX_opener) (PandaLogProxy logger, const char *args, PandaPreCheck *precheck, PandaCheck *check, void **user_data, PandaDestroy *destroy)
#        define HELP(desc, usage) const char *PANDACONCAT(PANDASEQ_MODULE,_LTX_desc) = desc; const char *PANDACONCAT(PANDASEQ_MODULE,_LTX_usage) = usage
#        define VER_INFO(version) const char *PANDACONCAT(PANDASEQ_MODULE,_LTX_version) = version
	EXTERN_C_END
#endif
