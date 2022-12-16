import { Router } from 'express';
import { errorCatcher, sanitize, methodNotAllowed, escapeRegExp, withRegex } from '../../helpers';
import { Database } from '../../Entities/CouchHelper';
import nano = require('nano');
import { User } from '../../Entities/entities';
import SearchWorker from '../../search_worker';
import { MAINTENANCE, MOLECULE_ROOT_DIR } from '../../constants';
import JSZip = require('jszip');
import fs from "fs";


const ListMoleculeRouter = Router();

interface Filters {
  /* Filters */
  version?: string;
  force_fields?: string;
  q: string;
  author: string;
  categories?: string;
  create_ways?: string;
  name?: string;
  alias?: string;
  owner?: string;

  /* Search modificators */
  is_regex?: string;
  combine?: "0" | "false" | "1" | "true";
  skip?: string;
  limit?: string;
  from_stashed?: "true" | "false"

}

ListMoleculeRouter.get('/', (req, res) => {
  (async () => {

    if (MAINTENANCE.mode) {
      res.json({ 'maintenance': true })
    } else {
      const query: nano.MangoQuery = { selector: {}, limit: 25, skip: 0 };

      const {
        force_fields,
        version,
        q: free_text,
        author,
        owner,
        categories,
        is_regex,
        create_ways,
        name: molecule_name,
        alias,
        combine: __as_version_tree__,
        skip,
        limit,
        from_stashed
      } = req.query as Partial<Filters>;

      console.log("from_stashed", from_stashed)
      const with_regex = is_regex === "true" || is_regex === "1";
      const bulk_request = __as_version_tree__ === "false" || __as_version_tree__ === "0";
      const search_in_stashed = from_stashed === "true"

      const selectors: any[] = [];
      if (force_fields) {
        const ff = force_fields.split(',');
        selectors.push({
          force_field: {
            $in: ff
          }
        });
      }
      if (create_ways) {
        const vers = create_ways.split(',');

        selectors.push({
          create_way: {
            $in: vers
          }
        });
      }
      if (version) {
        selectors.push({
          version: withRegex(version, with_regex)
        });
      }
      if (molecule_name) {
        selectors.push({
          name: withRegex(molecule_name, with_regex)
        });
      }
      if (alias) {
        selectors.push({
          $or: [
            { alias: withRegex(alias, with_regex) },
            {
              alternative_alias: {
                $elemMatch: withRegex(alias, with_regex)
              }
            }
          ]
        })
        // selectors.push({
        //   alias: withRegex(alias, with_regex)
        // });
        // selectors.push({
        //   alternative_alias : {
        //     $elemMatch : {
        //       $eq : withRegex(alias, with_regex)
        //     }
        //   }
        // })
      }
      if (author) {
        const users_that_match_query = await Database.user.find({
          limit: 99999,
          selector: {
            $or: [
              // @ts-ignore
              { email: withRegex(author, with_regex) },
              { name: withRegex(author, with_regex) },
            ]
          }
        }) as User[];

        const users_id = users_that_match_query.map(u => u.id);

        selectors.push({
          owner: { $in: users_id }
        });
      }
      if (owner) {
        selectors.push({
          owner
        });
      }
      if (categories) {
        const cat = categories.split(',');
        selectors.push({
          category: {
            $in: cat
          }
        });
      }
      if (free_text) {
        const search_text = with_regex ? author : escapeRegExp(free_text);
        const search_obj = {
          $regex: '(?i)' + search_text
        };

        selectors.push({
          $or: [
            { category: search_obj },
            { comments: search_obj },
            { name: search_obj },
            { alias: search_obj },
            { command_line: search_obj },
            { force_field: search_obj },
            { alternative_alias: { $elemMatch: search_obj } }
          ]
        });
      }

      if (selectors.length) {
        query.selector = { $and: selectors };
      }

      if (skip) {
        if (Number(skip) > 0) {
          query.skip = Number(skip);
        }
      }
      if (limit) {
        const l = Number(limit);

        if (l > 0 && l <= 200) {
          query.limit = l;
        }
      }


      const response = await SearchWorker.query(query, bulk_request);
      response.molecules = response.molecules.map(e => sanitize(e));


      // GET ITP FILE IN STRING
      if (response.length !== 0) {
        const file_id = response.molecules[0]["files"]

        const zipname = MOLECULE_ROOT_DIR + file_id + ".zip"
        fs.readFile(zipname, function (err, data) {
          if (err) throw err;
          JSZip.loadAsync(data).then(function (zip) {
            console.log(Object.keys(zip.files))
            const itpname = Object.keys(zip.files).filter((e) => e.endsWith(".itp"))[0]
            if (itpname !== null) {
              zip.file(itpname)!.async("string").then(function (data) {
                // data is a string
                // TODO your code goes here, this is just an example
                //console.log("blabla", data);
              })
            }
            else console.log(itpname, "NUUUL")

          });
        });
      }


      res.json(response);
    }

  })().catch(errorCatcher(res));
});

ListMoleculeRouter.all('/', methodNotAllowed('GET'))

export default ListMoleculeRouter;
